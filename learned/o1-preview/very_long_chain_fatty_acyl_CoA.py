"""
Classifies: CHEBI:61910 very long-chain fatty acyl-CoA
"""
"""
Classifies: CHEBI:51953 very long-chain fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_very_long_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a very long-chain fatty acyl-CoA based on its SMILES string.
    A very long-chain fatty acyl-CoA is a fatty acyl-CoA in which the fatty acyl group has a chain length greater than C22.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a very long-chain fatty acyl-CoA, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the Coenzyme A (CoA) moiety
    # Use a simplified SMARTS pattern for CoA identification (adenosine part)
    coa_smarts = '[C@@H]1O[C@H]([C@H](O)[C@@H]1OP(=O)(O)O)N2C=NC3=C2N=CN=C3N'
    coa_pattern = Chem.MolFromSmarts(coa_smarts)
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "CoA moiety not found"

    # Find the thioester linkage: C(=O)-S
    thioester_smarts = 'C(=O)S'
    thioester_pattern = Chem.MolFromSmarts(thioester_smarts)
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "Thioester linkage not found"

    # Assume there is only one acyl chain attached via thioester bond
    thioester_match = thioester_matches[0]
    carbonyl_c_idx = thioester_match[0]
    sulfur_idx = thioester_match[2]

    # Get the carbon atom attached to the carbonyl carbon that is not the oxygen or sulfur
    carbonyl_carbon = mol.GetAtomWithIdx(carbonyl_c_idx)
    acyl_chain_atom = None
    for neighbor in carbonyl_carbon.GetNeighbors():
        neighbor_idx = neighbor.GetIdx()
        if neighbor_idx != sulfur_idx and neighbor.GetAtomicNum() == 6:
            # This is the start of the acyl chain
            acyl_chain_atom = neighbor
            break

    if acyl_chain_atom is None:
        return False, "Acyl chain not found"

    # Traverse the acyl chain to count the number of carbons
    acyl_chain_carbon_idxs = set()
    visited_atoms = set()
    atoms_to_visit = [acyl_chain_atom.GetIdx()]
    while atoms_to_visit:
        atom_idx = atoms_to_visit.pop()
        if atom_idx in visited_atoms:
            continue
        visited_atoms.add(atom_idx)
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom.GetAtomicNum() == 6:  # Carbon atom
            acyl_chain_carbon_idxs.add(atom_idx)
            for neighbor in atom.GetNeighbors():
                neighbor_idx = neighbor.GetIdx()
                if neighbor_idx not in visited_atoms:
                    neighbor_atom = mol.GetAtomWithIdx(neighbor_idx)
                    if neighbor_atom.GetAtomicNum() == 6:
                        # Continue traversing carbons
                        atoms_to_visit.append(neighbor_idx)
        else:
            # Stop traversal at heteroatoms
            continue

    chain_length = len(acyl_chain_carbon_idxs)
    if chain_length > 22:
        return True, f"Acyl chain length is {chain_length}, which is greater than 22 carbons"
    else:
        return False, f"Acyl chain length is {chain_length}, which is not greater than 22 carbons"

__metadata__ = {   
    'chemical_class': {   
        'id': 'CHEBI:51953',
        'name': 'very long-chain fatty acyl-CoA',
        'definition': 'A fatty acyl-CoA in which the fatty acyl group has a chain length greater than C22.',
        'parents': ['CHEBI:37554']},
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
    'accuracy': None
}