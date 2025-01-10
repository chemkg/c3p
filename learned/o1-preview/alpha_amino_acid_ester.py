"""
Classifies: CHEBI:46874 alpha-amino acid ester
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_alpha_amino_acid_ester(smiles: str):
    """
    Determines if a molecule is an alpha-amino acid ester based on its SMILES string.
    An alpha-amino acid ester is an amino acid derivative obtained by esterification of the carboxyl group of an alpha-amino acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alpha-amino acid ester, False otherwise
        str: Reason for classification
    """
    # Parse SMILES to molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for alpha-amino acid ester
    # Pattern explanation:
    # [NX3;H2,H1;!$(NC=O)] - primary or secondary amine (not amide)
    # [CX4H] - tetrahedral carbon with one hydrogen (alpha carbon)
    # [C](=O)[O][#6] - esterified carboxyl group
    alpha_amino_acid_ester_smarts = "[NX3;H2,H1;!$(NC=O)][CX4H]([#6])[C](=O)[O][#6]"
    pattern = Chem.MolFromSmarts(alpha_amino_acid_ester_smarts)

    if pattern is None:
        return False, "Invalid SMARTS pattern"

    # Search for the alpha-amino acid ester pattern in the molecule
    matches = mol.GetSubstructMatches(pattern)
    if matches:
        return True, "Contains alpha-amino acid ester moiety"
    else:
        return False, "Does not contain alpha-amino acid ester moiety"

__metadata__ = {
    'chemical_class': {
        'id': None,
        'name': 'alpha-amino acid ester',
        'definition': 'The amino acid ester derivative obtained the formal condensation of an alpha-amino acid with an alcohol.',
        'parents': None
    },
    'config': {
        'llm_model_name': 'YourModelName',
        'f1_threshold': None,
        'max_attempts': None,
        'max_positive_instances': None,
        'max_positive_to_test': None,
        'max_negative_to_test': None,
        'max_positive_in_prompt': None,
        'max_negative_in_prompt': None,
        'max_instances_in_prompt': None,
        'test_proportion': None
    },
    'message': None,
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'precision': None,
    'recall': None,
    'f1': None,
    'accuracy': None
}