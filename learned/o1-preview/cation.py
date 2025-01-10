"""
Classifies: CHEBI:36916 cation
"""
"""
Classifies: CHEBI:36916 cation
"""
from rdkit import Chem

def is_cation(smiles: str):
    """
    Determines if a molecule is a cation based on its SMILES string.
    A cation is a monoatomic or polyatomic species having one or more elementary charges of the proton.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cation, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string into an RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Calculate the total formal charge of the molecule
    formal_charge = sum([atom.GetFormalCharge() for atom in mol.GetAtoms()])
    
    # Check if the molecule has a net positive charge
    if formal_charge > 0:
        return True, f"Molecule has net positive charge of {formal_charge}"
    else:
        return False, f"Molecule has net charge of {formal_charge}"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:36916',
        'name': 'cation',
        'definition': 'A monoatomic or polyatomic species having one or more elementary charges of the proton.',
        'parents': []
    },
    'config': {
        'llm_model_name': 'YourModelName',
        'f1_threshold': 0.8,
        'max_attempts': 5,
        'max_positive_instances': None,
        'max_positive_to_test': None,
        'max_negative_to_test': None,
        'max_positive_in_prompt': 50,
        'max_negative_in_prompt': 20,
        'max_instances_in_prompt': 100,
        'test_proportion': 0.1
    },
    'message': None,
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    # Additional metadata can be added here if necessary
}