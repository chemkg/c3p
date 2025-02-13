"""
Classifies: CHEBI:35692 dicarboxylic acid
"""
"""
Classifies: Dicarboxylic Acid – Any carboxylic acid containing two carboxy groups.
"""
from rdkit import Chem

def is_dicarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a dicarboxylic acid based on its SMILES string.
    A dicarboxylic acid is defined as any molecule containing exactly two carboxyl groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a dicarboxylic acid, False otherwise.
        str: Reason explaining the classification.
    """
    # Attempt to parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern for a carboxyl group.
    # The pattern "C(=O)[O;H1]" matches a carbon that is double-bonded to oxygen
    # and single-bonded to an oxygen carrying one hydrogen (i.e. an -OH group).
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)[O;H1]")
    if carboxyl_pattern is None:
        return False, "Error creating SMARTS for carboxyl group"
    
    # Get the number of carboxyl groups in the molecule.
    matches = mol.GetSubstructMatches(carboxyl_pattern)
    n_matches = len(matches)
    
    # Check if the molecule has exactly two carboxyl groups.
    if n_matches == 2:
        return True, "Molecule contains exactly two carboxyl groups; qualifies as a dicarboxylic acid."
    elif n_matches < 2:
        return False, f"Molecule contains only {n_matches} carboxyl group(s); need exactly 2."
    else:
        return False, f"Molecule contains {n_matches} carboxyl groups; exceeds the definition of a dicarboxylic acid."
        
# Example usage:
# To test the function, you can run:
# print(is_dicarboxylic_acid("OC(=O)C(O)=O"))  # Should return True for 2-formylglutaric acid or similar.