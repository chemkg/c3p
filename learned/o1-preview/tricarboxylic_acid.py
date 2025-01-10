"""
Classifies: CHEBI:27093 tricarboxylic acid
"""
"""
Classifies: CHEBI:39157 tricarboxylic acid
"""
from rdkit import Chem

def is_tricarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a tricarboxylic acid based on its SMILES string.
    A tricarboxylic acid is an oxoacid containing three carboxy groups (-COOH).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tricarboxylic acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define carboxylic acid SMARTS pattern
    carboxylic_acid_pattern = Chem.MolFromSmarts('[CX3](=O)[OX1H]')
    if carboxylic_acid_pattern is None:
        return False, "Invalid carboxylic acid SMARTS pattern"

    # Find carboxylic acid groups
    matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    num_carboxylic_acids = len(matches)

    # Check if there are exactly three carboxylic acid groups
    if num_carboxylic_acids == 3:
        return True, "Contains exactly three carboxylic acid groups"
    else:
        return False, f"Contains {num_carboxylic_acids} carboxylic acid group(s), expected exactly 3"

__metadata__ = {   'chemical_class': {   'id': 'CHEBI:39157',
                              'name': 'tricarboxylic acid',
                              'definition': 'An oxoacid containing three carboxy groups.',
                              'parents': ['CHEBI:33575']},
    'config': {   'llm_model_name': 'lbl/claude-sonnet',
                  'f1_threshold': 0.8,
                  'max_attempts': 5,
                  'max_positive_instances': None,
                  'max_positive_to_test': None,
                  'max_negative_to_test': None,
                  'max_positive_in_prompt': 50,
                  'max_negative_in_prompt': 20,
                  'max_instances_in_prompt': 100,
                  'test_proportion': 0.1},
    'message': None,
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': None,
    'num_false_positives': None,
    'num_true_negatives': None,
    'num_false_negatives': None,
    'num_negatives': None,
    'precision': None,
    'recall': None,
    'f1': None,
    'accuracy': None}