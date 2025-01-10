"""
Classifies: CHEBI:35819 branched-chain fatty acid
"""
"""
Classifies: CHEBI:59826 branched-chain fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_branched_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a branched-chain fatty acid based on its SMILES string.
    A branched-chain fatty acid is a fatty acid with one or more alkyl substituents on the parent hydrocarbon chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a branched-chain fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group (-C(=O)O)
    carboxylic_acid_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2H1]")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # Check for a long hydrocarbon chain (at least 6 carbons)
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 6:
        return False, "Carbon chain too short to be a fatty acid"

    # Check for alkyl substituents (branches)
    # A branch is defined as a carbon with at least 3 neighbors (sp3 hybridized)
    branch_pattern = Chem.MolFromSmarts("[CX4H3,CX4H2,CX4H1,CX4H0]")
    branch_matches = mol.GetSubstructMatches(branch_pattern)
    if len(branch_matches) < 2:  # At least one branch (excluding the carboxylic acid carbon)
        return False, "No alkyl substituents (branches) found"

    # Check if the branches are alkyl groups (methyl, ethyl, etc.)
    alkyl_branch_pattern = Chem.MolFromSmarts("[CX4H3,CX4H2,CX4H1,CX4H0]")
    alkyl_branch_matches = mol.GetSubstructMatches(alkyl_branch_pattern)
    if len(alkyl_branch_matches) < 2:  # At least one alkyl branch
        return False, "No alkyl substituents (branches) found"

    # Check for unsaturation (optional, as some branched-chain fatty acids can be unsaturated)
    # This is not strictly necessary but can be used to refine the classification
    double_bond_count = sum(1 for bond in mol.GetBonds() if bond.GetBondType() == Chem.BondType.DOUBLE)
    if double_bond_count > 0:
        return True, "Branched-chain fatty acid with unsaturation"
    else:
        return True, "Saturated branched-chain fatty acid"


__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:59826',
        'name': 'branched-chain fatty acid',
        'definition': 'Any fatty acid in which the parent hydrocarbon chain has one or more alkyl substituents; a common component in animal and bacterial lipids.',
        'parents': ['CHEBI:35366', 'CHEBI:26607']
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
    },
    'message': None,
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 150,
    'num_false_positives': 4,
    'num_true_negatives': 182407,
    'num_false_negatives': 23,
    'num_negatives': None,
    'precision': 0.974025974025974,
    'recall': 0.8670520231213873,
    'f1': 0.9174311926605504,
    'accuracy': 0.9998521228585199
}