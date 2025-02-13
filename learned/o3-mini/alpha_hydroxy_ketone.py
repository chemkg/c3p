"""
Classifies: CHEBI:139588 alpha-hydroxy ketone
"""
"""
Classifies: Alpha-hydroxy ketone
Definition: A ketone containing a hydroxy group on the alpha‐carbon relative to the C=O group.
Improved approach:
  • First, detect ketone groups (C=O with two carbon substituents).
  • For each ketone group, examine both alpha-carbons.
  • If an alpha-carbon is sp3-hybridized and has at least one oxygen substituent (which in turn bears one or more hydrogen atoms),
    then conclude the molecule contains an alpha‐hydroxy ketone.
This approach attempts to sidestep oversimplified SMARTS matches that can lead to false positives.
"""

from rdkit import Chem

def is_alpha_hydroxy_ketone(smiles: str):
    """
    Determines if a molecule is an alpha-hydroxy ketone based on its SMILES string.
    An alpha-hydroxy ketone is defined here as a ketone (C(=O) between two carbons) that has an -OH
    substituent on at least one of the carbons directly bonded to the carbonyl carbon.
    
    This improved method locates all ketone moieties (excluding aldehydes and acids)
    and then checks that one of the adjacent carbons (alpha positions) is sp3 and bears -OH.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as an alpha-hydroxy ketone.
        str: Explanation of the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern for a ketone: [#6]-[CX3](=O)-[#6]
    # This ensures the carbonyl carbon is attached to two carbon atoms.
    ketone_smarts = "[#6]-[CX3](=O)-[#6]"
    ketone_pattern = Chem.MolFromSmarts(ketone_smarts)
    ketone_matches = mol.GetSubstructMatches(ketone_pattern)
    
    if not ketone_matches:
        return False, "No ketone (C=O) group with two carbon substituents found."
    
    # For each ketone match, the pattern returns a tuple (idx_alpha, idx_C=O, idx_alpha2)
    # Check each of the two alpha-carbons for an attached hydroxy group.
    for match in ketone_matches:
        # The three indices for: alpha1, carbonyl, alpha2
        alpha1_idx, carbonyl_idx, alpha2_idx = match
        for alpha_idx in (alpha1_idx, alpha2_idx):
            alpha_atom = mol.GetAtomWithIdx(alpha_idx)
            # We require the alpha carbon to be sp3 (so that the OH is on a saturated carbon)
            if alpha_atom.GetHybridization() != Chem.rdchem.HybridizationType.SP3:
                continue  # skip if not sp3
            # Examine the neighbors of the alpha atom:
            for neighbor in alpha_atom.GetNeighbors():
                # Look for an oxygen atom.
                if neighbor.GetAtomicNum() == 8:
                    # Check if that oxygen has at least 1 hydrogen (it might be implicit though).
                    # RDKit’s GetTotalNumHs() returns total H count including implicit.
                    if neighbor.GetTotalNumHs() > 0:
                        # Found an OH group on an alpha carbon!
                        return True, "Molecule contains a ketone with an sp3 alpha-carbon bearing a hydroxy group (alpha-hydroxy ketone detected)."
    return False, "No ketone found with an sp3 alpha-carbon that bears a hydroxy substituent."
    
# Example usage (if desired):
if __name__ == "__main__":
    test_smiles = [
        "OCC(O)C(=O)CO",  # erythrulose, expected true
        "O=C(OCCC1=CC=C(O)C=C1)CCC(=O)[C@@H](O)C"  # Trichodenol B, expected false in the challenge
    ]
    for smi in test_smiles:
        res, reason = is_alpha_hydroxy_ketone(smi)
        print(f"SMILES: {smi}\nResult: {res}\nReason: {reason}\n")