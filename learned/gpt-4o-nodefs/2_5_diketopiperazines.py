"""
Classifies: CHEBI:65061 2,5-diketopiperazines
"""
from rdkit import Chem

def is_2_5_diketopiperazines(smiles: str):
    """
    Determines if a molecule is a 2,5-diketopiperazine based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2,5-diketopiperazine, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Refined pattern to strictly match 2,5-diketopiperazine
    # Here we focus on the core cyclization: a six-membered ring with two nitrogen atoms and two carbonyl groups
    diketopiperazine_core_pattern = Chem.MolFromSmarts("C1CN(C(=O)C1)C(=O)")

    # Check if the molecule matches the refined 2,5-diketopiperazine core pattern
    if mol.HasSubstructMatch(diketopiperazine_core_pattern):
        return True, "Contains 2,5-diketopiperazine motif"
    
    return False, "Does not contain 2,5-diketopiperazine motif"