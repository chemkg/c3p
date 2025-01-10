"""
Classifies: CHEBI:24128 furanocoumarin
"""
"""
Classifies: furanocoumarin
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_furanocoumarin(smiles: str):
    """
    Determines if a molecule is a furanocoumarin based on its SMILES string.
    A furanocoumarin is a furochromene that consists of a furan ring fused with a coumarin.
    The fusion may occur in different ways to give several isomers.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a furanocoumarin, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # General SMARTS pattern for furanocoumarin core (tricyclic ring system)
    furanocoumarin_pattern = Chem.MolFromSmarts('c1cc2c(c1)oc(=O)c3ccoc23')

    if mol.HasSubstructMatch(furanocoumarin_pattern):
        return True, "Contains furanocoumarin core"

    else:
        # Check for angular furanocoumarin core (angelicin-type)
        angular_pattern = Chem.MolFromSmarts('c1ccc2c(c1)oc(=O)c3ccoc23')
        if mol.HasSubstructMatch(angular_pattern):
            return True, "Contains angular furanocoumarin core"
        else:
            return False, "Does not contain furanocoumarin core"

__metadata__ = {
    'chemical_class': {
        'id': None,
        'name': 'furanocoumarin',
        'definition': 'Any furochromene that consists of a furan ring fused with a coumarin. The fusion may occur in different ways to give several isomers.',
        'parents': []
    },
    'config': {
        'llm_model_name': 'lbl/claude-sonnet',
        'f1_threshold': 0.8,
        'max_attempts': 5,
        'max_positive_instances': None,
        'max_positive_to_test': None,
        'max_negative_to_test': None,
        'max_positive_in_prompt': 50,
        'max_negative_in_prompt': 20,
        'max_instances_in_prompt': 100,
        'test_proportion': 0.1
    }
}