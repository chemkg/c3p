"""
Classifies: CHEBI:75659 O-acyl-L-carnitine
"""
"""
Classifies: O-acyl-L-carnitine (an O-acylcarnitine in which the carnitine component has L-configuration)

We use the following approach:
  1. Parse the SMILES.
  2. Look for a substructure matching the carnitine moiety with an O-acyl substituent.
     In our examples, this substructure is roughly:
       a chiral carbon (mapped as atom 1) attached to:
         - an oxygen that in turn is bonded to a carbonyl (the acyl group),
         - a propanoate (carboxylate) chain,
         - and a trimethylammonium moiety.
  3. Then, after assigning stereochemistry using RDKit, we check that the CIP label
     at the mapped chiral center is "R" (which corresponds to L-carnitine).
     
If the molecule does not contain the expected pattern, the function returns False
and an explanation.
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_O_acyl_L_carnitine(smiles: str):
    """
    Determines if a molecule is an O-acyl-L-carnitine based on its SMILES string.
    An O-acyl-L-carnitine must contain a carnitine skeleton with an acyl group esterified
    to its hydroxyl group, and the carnitine chiral center must have the natural L-configuration,
    which in CIP nomenclature for carnitine is “R”.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is an O-acyl-L-carnitine, False otherwise.
        str: Reason for the classification decision.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Assign stereochemistry so that the CIP code (_CIPCode property) can be computed
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    
    # Define two SMARTS patterns to cover both possible chirality annotations [C@H] and [C@@H]
    # The pattern seeks a chiral carbon (mapped as atom 1) that is bound to:
    #   - [O][C](=O)[*]   : an oxygen connected to a carbonyl (the acyl substituent; [*] allows any acyl chain)
    #   - C[N+](C)(C)C   : a trimethylammonium group
    #   - CC(=O)[O-]     : a carboxylate portion (the propionate side-chain)
    pattern_smarts1 = "[C@H:1]([O][C](=O)[*])(C[N+](C)(C)C)(CC(=O)[O-])"
    pattern_smarts2 = "[C@@H:1]([O][C](=O)[*])(C[N+](C)(C)C)(CC(=O)[O-])"
    
    patt1 = Chem.MolFromSmarts(pattern_smarts1)
    patt2 = Chem.MolFromSmarts(pattern_smarts2)
    
    matches = []
    if patt1:
        matches.extend(mol.GetSubstructMatches(patt1, useChirality=True))
    if patt2:
        matches.extend(mol.GetSubstructMatches(patt2, useChirality=True))
    
    if not matches:
        return False, "No substructure matching an O-acyl-carnitine skeleton was found"
    
    # Loop over all matches. For each, check the stereochemistry at the mapped chiral atom.
    # Our SMARTS was defined so that the first atom in the query (index 0 in the match tuple) is the chiral carbon.
    for match in matches:
        # Get the chiral center atom (from the molecule) corresponding to the mapped atom (map number 1)
        chiral_atom = mol.GetAtomWithIdx(match[0])
        # Ensure that the CIP code property is present for this atom.
        if not chiral_atom.HasProp("_CIPCode"):
            continue  # try next match if CIP info is missing
        cip = chiral_atom.GetProp("_CIPCode")
        # In carnitine, the naturally occurring L-form has CIP code "R". (L-carnitine is (R)-carnitine.)
        if cip == "R":
            return True, "Matches O-acyl-L-carnitine pattern with correct (R) stereochemistry at carnitine center"
        else:
            # The acylcarnitine motif is present but the chiral center is not the L (R) configuration.
            return False, f"carnitine center has CIP code {cip} (expected 'R' for L-carnitine)"
    
    return False, "No match found with defined stereochemistry information"

# Example usage:
if __name__ == '__main__':
    # Test with one of the provided examples: O-acetyl-L-carnitine
    test_smiles = "CC(=O)O[C@H](CC([O-])=O)C[N+](C)(C)C"
    result, reason = is_O_acyl_L_carnitine(test_smiles)
    print("Result:", result)
    print("Reason:", reason)