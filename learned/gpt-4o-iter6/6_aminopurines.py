"""
Classifies: CHEBI:20706 6-aminopurines
"""
from rdkit import Chem

def is_6_aminopurines(smiles: str):
    """
    Determines if a molecule belongs to the 6-aminopurines class based on its SMILES string.
    6-aminopurine (adenine) is characterized by a specific purine structure with an amino group at the 6th position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule contains 6-aminopurine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define adenine substructure patterns considering common variations of adenine
    adenine_patterns = [
        Chem.MolFromSmarts("n1cnc2c(ncn2)n1"),  # Basic adenine
        Chem.MolFromSmarts("c1[nH]c2[nH]cnc2n1"),  # Alternative tautomer
        Chem.MolFromSmarts("n1c[nH]c2c1ncn2")    # Another common form
    ]
    
    # Search for adenine substructure in the molecule
    for pattern in adenine_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains 6-aminopurine (adenine) as part of its structure"
    
    return False, "Does not contain 6-aminopurine (adenine) structure"