"""
Classifies: CHEBI:139588 alpha-hydroxy ketone
"""
"""
Classifies: Alpha-hydroxy ketone
Definition: A ketone containing a hydroxy group on the alpha‐carbon relative to the C=O group.
Improved approach:
  • Locate ketone moieties using a SMARTS pattern ensuring the carbonyl carbon is sandwiched between two carbons.
  • For each ketone match, ensure the match tuple exactly contains three atoms.
  • Check if one of the alpha-carbons (first and third atoms) is sp3 and has an oxygen neighbor bearing at least one hydrogen.
"""

from rdkit import Chem

def is_alpha_hydroxy_ketone(smiles: str):
    """
    Determines if a molecule is an alpha-hydroxy ketone based on its SMILES string.
    An alpha-hydroxy ketone is defined as a ketone (C(=O) linked to two carbons) with an -OH 
    substituent on at least one of the carbons adjacent to the C=O group.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is an alpha-hydroxy ketone, False otherwise.
        str: Explanation for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern for a ketone with two carbon substituents:
    #    [#6]-[CX3](=O)-[#6]
    # This ensures that the carbonyl carbon (with =O) is flanked by two carbon atoms.
    ketone_smarts = "[#6]-[CX3](=O)-[#6]"
    ketone_pattern = Chem.MolFromSmarts(ketone_smarts)
    ketone_matches = mol.GetSubstructMatches(ketone_pattern)
    
    if not ketone_matches:
        return False, "No ketone (C=O) group with two carbon substituents found."
    
    # Iterate over each ketone match.
    # We expect each match to be a tuple of 3 atom indices: (alpha1, carbonyl, alpha2).
    for match in ketone_matches:
        # If the match does not contain exactly three atoms, skip it.
        if len(match) != 3:
            continue
        alpha1_idx, carbonyl_idx, alpha2_idx = match
        
        # Check both alpha-carbons for a hydroxy group.
        for alpha_idx in (alpha1_idx, alpha2_idx):
            alpha_atom = mol.GetAtomWithIdx(alpha_idx)
            # Require the alpha carbon to be sp3 hybridized.
            if alpha_atom.GetHybridization() != Chem.rdchem.HybridizationType.SP3:
                continue
            # Now, check the neighbors of the alpha carbon for an oxygen atom.
            for neighbor in alpha_atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 8:
                    # Check if the oxygen bears at least one hydrogen (using total H count, including implicit).
                    if neighbor.GetTotalNumHs() > 0:
                        return True, ("Molecule contains a ketone with an sp3 alpha-carbon bearing a hydroxy group "
                                      "(alpha-hydroxy ketone detected).")
    
    return False, "No ketone found with an sp3 alpha-carbon that bears a hydroxy substituent."

# Example usage:
if __name__ == "__main__":
    # Test examples; you may extend these or replace with your own test cases.
    test_smiles = [
        "OCC(O)C(=O)CO",  # erythrulose, should classify as alpha-hydroxy ketone (expected True)
        "C[C@H]1[C@H]2[C@H](Cc3c[nH]c4ccccc34)NC(=O)[C@@]22[C@@H](\\C=C\\C[C@H](C)\\C=C(C)\\C(=O)[C@@H](O)CCC2=O)[C@@H]2O[C@]12C"  # chaetoglobosin F, complex structure, testing error avoidance.
    ]
    for smi in test_smiles:
        res, reason = is_alpha_hydroxy_ketone(smi)
        print(f"SMILES: {smi}\nResult: {res}\nReason: {reason}\n")