"""
Classifies: CHEBI:39410 1,2,4-triazines
"""
from rdkit import Chem

def is_1_2_4_triazines(smiles: str):
    """
    Determines if a molecule is a 1,2,4-triazine based on its SMILES string.
    A 1,2,4-triazine is a six-membered heterocyclic ring containing three nitrogens at
    positions 1, 2, and 4.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 1,2,4-triazine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for 1,2,4-triazine ring pattern
    triazine_pattern = Chem.MolFromSmarts("n1cnncn1")
    if mol.HasSubstructMatch(triazine_pattern):
        return True, "Contains a 1,2,4-triazine ring"

    return False, "Does not contain a 1,2,4-triazine ring"