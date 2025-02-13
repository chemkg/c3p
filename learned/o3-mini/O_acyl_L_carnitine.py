"""
Classifies: CHEBI:75659 O-acyl-L-carnitine
"""
"""
Classifies: O-acyl-L-carnitine
Definition: An O-acylcarnitine in which the carnitine component has L-configuration.
The molecule is expected to contain an ester linkage between an acyl group and a carnitine backbone.
The carnitine substructure is identified by a chiral carbon (mapped as :1) 
attached to an ester oxygen (mapped as :2), the acyl carbon (mapped as :3),
a trimethylammonium fragment (mapped as :4) and a propanoate (carboxylate) group.
In addition, the CIP stereochemistry computed on the carnitine chiral center must be "R" (which is typical for L–carnitine).
We then do a DFS on the acyl portion (starting at the acyl carbon but not crossing back via the bridging oxygen)
to check that if the acyl fragment is very short (<6 heavy atoms) it does not have an extra carboxylate group 
(which would be unusual).
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_O_acyl_L_carnitine(smiles: str):
    """
    Determine if a molecule is an O-acyl-L-carnitine based on its SMILES string.
    It must contain an O-acyl–carnitine skeleton in which the carnitine chiral
    center (mapped as :1) has CIP code "R". We also attempt to reject cases where
    the acyl substituent (attached via the oxygen, mapped as :2 and :3) shows an extra carboxylate 
    group for a short chain (which would be atypical).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule meets the criteria, False otherwise.
        str: Explanation for the decision.
    """
    # Parse the SMILES string to a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Compute stereochemistry so that CIP codes are assigned to chiral centers.
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    
    # Define two SMARTS patterns to capture the two possibilities of chirality ([C@H] and [C@@H])
    # The SMARTS looks for:
    #    - a chiral carbon (mapping :1)
    #    - bearing an oxygen attached to an acyl carbon (mapping :2 and :3 respectively)
    #    - a trimethylammonium fragment (mapping :4)
    #    - and a propanoate branch (non-mapped).
    pattern_smarts1 = "[C@H:1]([O:2][C:3](=O))([C:4][N+](C)(C)C)(CC(=O)[O-])"
    pattern_smarts2 = "[C@@H:1]([O:2][C:3](=O))([C:4][N+](C)(C)C)(CC(=O)[O-])"
    patt1 = Chem.MolFromSmarts(pattern_smarts1)
    patt2 = Chem.MolFromSmarts(pattern_smarts2)
    if patt1 is None or patt2 is None:
        return False, "Could not compile SMARTS patterns"
    
    # Gather all matches using chirality in matching
    matches = []
    matches.extend(mol.GetSubstructMatches(patt1, useChirality=True))
    matches.extend(mol.GetSubstructMatches(patt2, useChirality=True))
    
    if not matches:
        return False, "No substructure matching an O-acyl–carnitine skeleton was found"
    
    # A helper function to perform a depth-first-search (DFS) on the molecule graph
    # starting from 'start_idx' and avoiding a specific atom (exclude_idx, the bridging oxygen).
    def dfs_fragment(mol, start_idx, exclude_idx):
        visited = set()
        stack = [start_idx]
        while stack:
            curr = stack.pop()
            if curr in visited:
                continue
            visited.add(curr)
            atom = mol.GetAtomWithIdx(curr)
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() == exclude_idx:
                    continue  # do not cross the bond to the bridging oxygen
                if nbr.GetIdx() not in visited:
                    stack.append(nbr.GetIdx())
        return visited

    # Define a SMARTS for a deprotonated carboxylate group.
    ca_pattern = Chem.MolFromSmarts("[CX3](=O)[O-]")
    
    valid_match_found = False
    for match in matches:
        # According to our mapping:
        # match[0] -> carnitine chiral center
        # match[1] -> bridging oxygen (ester oxygen)
        # match[2] -> acyl carbon (carbonyl carbon of the acyl group)
        # match[3] -> carbon adjacent to the trimethylammonium
        try:
            carnitine_center = mol.GetAtomWithIdx(match[0])
            bridging_oxygen = mol.GetAtomWithIdx(match[1])
            acyl_carbon = mol.GetAtomWithIdx(match[2])
        except IndexError:
            continue
        
        # Check that the chiral center has CIP code "R"
        if not carnitine_center.HasProp("_CIPCode"):
            continue
        if carnitine_center.GetProp("_CIPCode") != "R":
            continue  # wrong stereochemistry for L-carnitine
        
        # Now check the acyl portion.
        # We do a DFS starting at the acyl carbon and disallow going back via the bridging oxygen.
        acyl_fragment_indices = dfs_fragment(mol, match[2], match[1])
        # Count heavy atoms (non-hydrogen atoms) in this acyl fragment.
        heavy_atoms = [idx for idx in acyl_fragment_indices if mol.GetAtomWithIdx(idx).GetAtomicNum() > 1]
        
        # The idea is that if the acyl fragment contains an extra carboxylate (beyond what is
        # expected in a simple acyl chain) and its heavy atom count is small (< 6), then the acyl group
        # does not represent a typical fatty-acyl moiety.
        # We create a SMILES of the fragment and then a new Mol from it.
        frag_smiles = Chem.MolFragmentToSmiles(mol, atomsToUse=list(acyl_fragment_indices), canonical=True)
        acyl_frag = Chem.MolFromSmiles(frag_smiles) if frag_smiles else None
        
        extra_check_pass = True
        if acyl_frag is not None:
            # Look for carboxylate pattern(s) in the acyl fragment.
            # Note that in a regular acyl chain the ester carbonyl oxygen is not charged,
            # so the presence of a deprotonated carboxylate ([O-]) suggests a diacid.
            ca_matches = acyl_frag.GetSubstructMatches(ca_pattern)
            if ca_matches and len(heavy_atoms) < 6:
                extra_check_pass = False
        
        # If the acyl fragment passes our extra check, then we have found a valid O-acyl-L-carnitine.
        if extra_check_pass:
            valid_match_found = True
            break

    if valid_match_found:
        return True, "Matches O-acyl-L-carnitine pattern with correct (R) stereochemistry at the carnitine center"
    else:
        return False, "Found a carnitine-like fragment but either the CIP configuration was not 'R' or the acyl group was inconsistent with typical properties"

# A quick test section (this part will run if the file is executed as a script)
if __name__ == '__main__':
    # Example: O-acetyl-L-carnitine
    test_smiles = "CC(=O)O[C@H](CC([O-])=O)C[N+](C)(C)C"
    result, reason = is_O_acyl_L_carnitine(test_smiles)
    print("Result:", result)
    print("Reason:", reason)