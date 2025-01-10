"""
Classifies: CHEBI:27283 very long-chain fatty acid
"""
"""
Classifies: CHEBI:27388 very long-chain fatty acid
"""
from rdkit import Chem

def is_very_long_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a very long-chain fatty acid based on its SMILES string.
    A very long-chain fatty acid is defined as a fatty acid with a chain length greater than C22.
    Ultra-long-chain fatty acids have a chain length greater than C27.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a very long-chain fatty acid, False otherwise
        str: Reason for classification
    """
    from collections import deque

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for carboxylic acid group (-COOH or -COO-)
    carboxylic_acid = Chem.MolFromSmarts('[CX3](=O)[O;H1,-1]')
    matches = mol.GetSubstructMatches(carboxylic_acid)
    if len(matches) == 0:
        return False, "No carboxylic acid group found"
    elif len(matches) > 1:
        return False, "Multiple carboxylic acid groups found"

    # Get the carboxyl carbon atom index
    carboxyl_carbon_idx = matches[0][0]
    carboxyl_carbon = mol.GetAtomWithIdx(carboxyl_carbon_idx)

    # Initialize variables for BFS traversal
    max_chain_length = 0
    visited = set()

    # Function to perform BFS and find the longest carbon chain
    def bfs(start_atom_idx):
        visited = set()
        queue = deque()
        queue.append((start_atom_idx, 1))
        local_max = 1
        while queue:
            current_idx, length = queue.popleft()
            visited.add(current_idx)
            atom = mol.GetAtomWithIdx(current_idx)
            for neighbor in atom.GetNeighbors():
                neighbor_idx = neighbor.GetIdx()
                if neighbor_idx not in visited and neighbor.GetAtomicNum() == 6:
                    queue.append((neighbor_idx, length + 1))
                    if length + 1 > local_max:
                        local_max = length + 1
        return local_max

    # Start BFS traversal from the carboxyl carbon
    max_chain_length = bfs(carboxyl_carbon_idx)

    # Check if chain length exceeds 22 carbons
    if max_chain_length > 22:
        if max_chain_length > 27:
            return True, f"Ultra-long-chain fatty acid with chain length C{max_chain_length}"
        else:
            return True, f"Very long-chain fatty acid with chain length C{max_chain_length}"
    else:
        return False, f"Chain length is C{max_chain_length}, not greater than C22"

__metadata__ = {   'chemical_class': {   'id': 'CHEBI:27388',
                          'name': 'very long-chain fatty acid',
                          'definition': 'A fatty acid which has a chain length greater than C22. Very long-chain fatty acids which have a chain length greater than C27 are also known as ultra-long-chain fatty acids.',
                          'parents': []},
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