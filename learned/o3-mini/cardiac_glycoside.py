"""
Classifies: CHEBI:83970 cardiac glycoside
"""
"""
Classifies: Cardiac glycoside (Steroid lactones containing sugar residues that act on the contractile force of the cardiac muscles)
This heuristic approach searches for:
  1. A basic steroid nucleus (four fused rings typical of steroids).
  2. A lactone moiety (here a furanone ring pattern).
  3. At least one sugar residue (pyranose or furanose ring motif).

Note: The SMARTS patterns used are approximate and may not capture every variant.
"""

from rdkit import Chem

def is_cardiac_glycoside(smiles: str):
    """
    Determines if a molecule is a cardiac glycoside based on its SMILES string.
    Cardiac glycosides are defined as steroid lactones containing sugar residues.
    
    Heuristic criteria:
      - Steroid nucleus: a tetracyclic fused ring system.
      - Lactone ring: a five‐membered ring containing a carbonyl and an oxygen (butenolide).
      - Sugar residue: at least one sugar ring (pyranose or furanose pattern).
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule qualifies as a cardiac glycoside, False otherwise.
        str: Reason for the classification.
    """
    
    # Parse the input SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # --- Check for the steroid nucleus ---
    # A rough SMARTS for a steroid backbone (four fused rings) is used.
    # This pattern is only approximate.
    steroid_smarts = "C1CC2CCC3C(C2)CCC1C3"
    steroid_pattern = Chem.MolFromSmarts(steroid_smarts)
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid nucleus (four fused rings) detected"
    
    # --- Check for lactone ring ---
    # Many cardiac glycosides have a butenolide (furanone) moiety.
    # We look for a five-membered ring with a carbonyl (C=O) and an adjacent ring oxygen.
    lactone_smarts = "C1=CC(=O)O1"
    lactone_pattern = Chem.MolFromSmarts(lactone_smarts)
    if not mol.HasSubstructMatch(lactone_pattern):
        return False, "No lactone (butenolide) ring detected"
    
    # --- Check for sugar residues ---
    # We attempt to detect at least one sugar ring.
    # Two simple SMARTS patterns are used:
    #   a pyranose-like pattern: a six-membered ring with one ring oxygen and five carbons bearing hydroxyls.
    #   a furanose-like pattern: a five-membered ring with one ring oxygen and four carbons bearing hydroxyls.
    pyranose_smarts = "OC1C(O)C(O)C(O)C(O)C1"   # simplified (no stereochemistry)
    furanose_smarts = "OC1C(O)C(O)C(O)O1"
    pyranose_pattern = Chem.MolFromSmarts(pyranose_smarts)
    furanose_pattern = Chem.MolFromSmarts(furanose_smarts)
    
    sugar_found = False
    if mol.HasSubstructMatch(pyranose_pattern):
        sugar_found = True
    elif mol.HasSubstructMatch(furanose_pattern):
        sugar_found = True
        
    if not sugar_found:
        return False, "No sugar residue (glycosidic moiety) detected"
    
    # If all criteria are met, classify as cardiac glycoside.
    msg = ("Detected a steroid nucleus, a lactone (butenolide) ring, and sugar residue(s) "
           "consistent with a cardiac glycoside structure.")
    return True, msg

# Example usage:
if __name__ == "__main__":
    # Test with one example SMILES from the provided list (e.g., Erychroside)
    test_smiles = "O[C@@]12[C@]3([C@@]([C@@]4([C@](O)(CC3)C[C@@H](O[C@@H]5O[C@@H]([C@@H](O[C@@H]6OC[C@@H](O)[C@H](O)[C@H]6O)[C@@H](O)C5)C)CC4)C=O)(CC[C@@]1([C@H](CC2)C=7COC(=O)C7)C)[H])[H]"
    result, reason = is_cardiac_glycoside(test_smiles)
    print("Is cardiac glycoside?", result)
    print("Reason:", reason)