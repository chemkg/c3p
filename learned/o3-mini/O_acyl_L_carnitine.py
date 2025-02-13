"""
Classifies: CHEBI:75659 O-acyl-L-carnitine
"""
"""
Classifies: O-acyl-L-carnitine
Definition: An O-acylcarnitine in which the carnitine component has L-configuration.
We require that the molecule contains a carnitine-like fragment defined roughly as:
  - A chiral carbon (mapped as :1) attached to:
      • An oxygen atom that is esterified with an acyl group. In our match we include only up to 
        the acyl carbonyl carbon (mapped as :3).
      • A trimethylammonium group (mapped as :4).
      • A propanoate (carboxylate) group.
  - When CIP stereochemistry is computed the chiral center must have CIP code "R"
    (the natural configuration for L–carnitine).
In addition we try to filter out cases where the acyl part is “unusual” (e.g. extra carboxylate on a short chain).
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_O_acyl_L_carnitine(smiles: str):
    """
    Determine if a molecule is an O-acyl-L-carnitine based on its SMILES string.
    It must contain an O-acyl–carnitine skeleton in which the carnitine chiral
    center (the atom mapped as :1) has CIP code "R". We also attempt to reject cases where
    the acyl substituent (attached via the oxygen, mapped as :2 and :3) shows a short diacid behavior.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule meets the criteria, False otherwise.
        str: Explanation for the decision.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Compute stereochemistry so that CIP codes are assigned
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    
    # Define two SMARTS patterns covering the two possible chirality notations.
    # The pattern covers:
    #   - A chiral carbon (mapping :1) attached to
    #       • An oxygen (mapping :2) that is bonded to a carbonyl carbon (mapping :3) via a single bond;
    #         the carbonyl oxygen (as double-bond partner) is not mapped.
    #       • A trimethylammonium fragment [C:4][N+](C)(C)C
    #       • A propanoate branch written as CC(=O)[O-] (note that this branch is not mapped)
    #
    # By stopping the acyl part at the carbonyl (mapped as :3) we avoid incorporating further substructure.
    pattern_smarts1 = "[C@H:1]([O:2][C:3](=O))([C:4][N+](C)(C)C)(CC(=O)[O-])"
    pattern_smarts2 = "[C@@H:1]([O:2][C:3](=O))([C:4][N+](C)(C)C)(CC(=O)[O-])"
    patt1 = Chem.MolFromSmarts(pattern_smarts1)
    patt2 = Chem.MolFromSmarts(pattern_smarts2)
    
    if patt1 is None or patt2 is None:
        return False, "Could not compile SMARTS patterns"
    
    # Get all matches (using chirality in matching)
    matches = []
    matches.extend(mol.GetSubstructMatches(patt1, useChirality=True))
    matches.extend(mol.GetSubstructMatches(patt2, useChirality=True))
    
    if not matches:
        return False, "No substructure matching an O-acyl–carnitine skeleton was found"
    
    valid_match_found = False
    for match in matches:
        # match tuple ordering is by mapping number.
        # match[0] is the carnitine chiral center.
        try:
            carnitine_center = mol.GetAtomWithIdx(match[0])
        except IndexError:
            continue

        # Check if CIP code was assigned to the chiral center:
        if not carnitine_center.HasProp("_CIPCode"):
            continue
        cip = carnitine_center.GetProp("_CIPCode")
        if cip != "R":
            continue  # wrong configuration for L-carnitine
        
        # At this point the match has the expected stereochemistry.
        # In our mapping, match[1] is the bridging oxygen and match[2] is the acyl carbon.
        acyl_atom = mol.GetAtomWithIdx(match[2])
        # Let us see if the acyl chain “extends” beyond the carbonyl carbon.
        # (In many cases we expect the acyl part to be a fatty acid chain; however, if the acyl
        # portion contains an extra carboxylate coming off a short chain (e.g. succinyl), we wish to reject it.)
        ext_neighbors = []
        for nbr in acyl_atom.GetNeighbors():
            # Exclude the bridging oxygen (mapped as match[1])
            if nbr.GetIdx() == match[1]:
                continue
            # Exclude the double-bonded oxygen (the carbonyl oxygen) – its bond is DOUBLE.
            bond = mol.GetBondBetweenAtoms(acyl_atom.GetIdx(), nbr.GetIdx())
            if bond and bond.GetBondType() == Chem.BondType.DOUBLE and nbr.GetAtomicNum() == 8:
                continue
            ext_neighbors.append(nbr)
        # If there is at least one neighbor beyond those required for a simple ester,
        # then we try to “extract” the acyl fragment.
        extra_check_pass = True
        if ext_neighbors:
            # Make a copy so that we can break the ester bond.
            rmol = Chem.RWMol(mol)
            try:
                rmol.RemoveBond(match[1], match[2])
            except Exception:
                pass
            frags = Chem.GetMolFrags(rmol, asMols=True, sanitizeFrags=True)
            # Find the fragment containing the acyl carbon (match[2])
            acyl_frag = None
            for frag in frags:
                if frag.HasAtomWithIdx(match[2]):
                    acyl_frag = frag
                    break
            if acyl_frag is not None:
                # Count heavy atoms (exclude H)
                heavy_atoms = [a for a in acyl_frag.GetAtoms() if a.GetAtomicNum() > 1]
                # Look for an extra carboxylate group in this fragment (a deprotonated carbonyl)
                ca_pattern = Chem.MolFromSmarts("[CX3](=O)[O-]")
                extra_acid = acyl_frag.HasSubstructMatch(ca_pattern)
                # The idea is that if a carboxylate appears in the acyl fragment
                # and the fragment is short (say, fewer than 6 heavy atoms),
                # then it’s likely an undesired dicarboxylic acid (like succinyl).
                if extra_acid and len(heavy_atoms) < 6:
                    extra_check_pass = False
        if not extra_check_pass:
            continue
        
        # Passed all tests (correct CIP and no “short” diacid as acyl group)
        valid_match_found = True
        break

    if valid_match_found:
        return True, "Matches O-acyl-L-carnitine pattern with correct (R) stereochemistry at carnitine center"
    else:
        # If we found some matches but none had the correct stereochemistry or passed the acyl check,
        # we return an appropriate explanation.
        return False, "Found a carnitine-like fragment but either the CIP configuration was not 'R' or the acyl group was inconsistent with O–acyl–L–carnitine"

# Example usage (for quick test):
if __name__ == '__main__':
    # Try one example: O-acetyl-L-carnitine
    test_smiles = "CC(=O)O[C@H](CC([O-])=O)C[N+](C)(C)C"
    result, reason = is_O_acyl_L_carnitine(test_smiles)
    print("Result:", result)
    print("Reason:", reason)