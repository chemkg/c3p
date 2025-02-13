"""
Classifies: CHEBI:20706 6-aminopurines
"""
"""
Classifies: 6-aminopurines (any compound having 6-aminopurine (adenine) as part of its structure)
"""
from rdkit import Chem

def is_6_aminopurines(smiles: str):
    """
    Determines if a molecule contains a 6-aminopurine (adenine) moiety.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule contains the 6-aminopurine substructure, False otherwise.
        str: Detailed reason for the classification.
    """
    # Parse the input SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the adenine substructure.
    # Adenine is represented as 6-aminopurine, and
    # its SMILES is given by "Nc1ncnc2ncnc12".
    adenine = Chem.MolFromSmiles("Nc1ncnc2ncnc12")
    if adenine is None:
        return False, "Error creating adenine substructure query"
    
    # Check if the molecule contains the adenine substructure
    if mol.HasSubstructMatch(adenine):
        return True, "Molecule contains the 6-aminopurine (adenine) moiety"
    else:
        return False, "6-aminopurine (adenine) moiety not found in the molecule"
        
# Example usage (uncomment to test)
# example_smiles = "Cn1cnc(N)c2ncnc12"  # 3-methyladenine is an example that should return True.
# result, reason = is_6_aminopurines(example_smiles)
# print(result, reason)