"""
Classifies: CHEBI:28863 flavanones
"""
"""
Classifies: CHEBI:50502 flavanones
"""
from rdkit import Chem

def is_flavanones(smiles: str):
    """
    Determines if a molecule is a flavanone based on its SMILES string.
    A flavanone is characterized by the 2,3-dihydro-2-arylchromen-4-one core structure.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a flavanone, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define flavanone core SMARTS pattern (without stereochemistry)
    flavanone_smarts = 'O=C1CCC(Oc2ccccc2)C1'  # Core flavanone structure
    flavanone_pattern = Chem.MolFromSmarts(flavanone_smarts)
    if flavanone_pattern is None:
        return False, "Invalid flavanone core SMARTS pattern"

    # Check for flavanone core structure
    if not mol.HasSubstructMatch(flavanone_pattern):
        return False, "Flavanone core structure not found"

    return True, "Contains flavanone core structure (2,3-dihydro-2-arylchromen-4-one)"

__metadata__ = {   'chemical_class': {   'id': 'CHEBI:50502',
                                  'name': 'flavanones',
                                  'definition': 'Members of the class of flavans with a 3,4-dihydro-2-aryl-2H-1-benzopyran-4-one skeleton and its substituted derivatives.',
                                  'parents': ['CHEBI:135930', 'CHEBI:2455']},
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
            'attempt': 1,
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