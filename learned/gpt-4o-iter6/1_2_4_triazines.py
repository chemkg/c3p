"""
Classifies: CHEBI:39410 1,2,4-triazines
"""
from rdkit import Chem

def is_1_2_4_triazines(smiles: str):
    """
    Determines if a molecule is a 1,2,4-triazine based on its SMILES string.
    A 1,2,4-triazine is characterized by a six-membered aromatic ring with nitrogen atoms at positions 1, 2, and 4.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a 1,2,4-triazine, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the SMARTS pattern for the 1,2,4-triazine ring
    triazine_pattern = Chem.MolFromSmarts("n1cncnc1") # 1,2,4-triazine core SMARTS
    
    # Check if the molecule contains the 1,2,4-triazine substructure
    if mol.HasSubstructMatch(triazine_pattern):
        return True, "Contains the 1,2,4-triazine core structure"
    
    return False, "Does not contain the 1,2,4-triazine core structure"