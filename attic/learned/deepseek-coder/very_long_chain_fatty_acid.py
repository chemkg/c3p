"""
Classifies: CHEBI:27283 very long-chain fatty acid
"""
"""
Classifies: CHEBI:27283 very long-chain fatty acid
"""
from rdkit import Chem

def is_very_long_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a very long-chain fatty acid based on its SMILES string.
    A very long-chain fatty acid is defined as having a carbon chain length greater than C22.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a very long-chain fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group (-COOH)
    carboxylic_acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # Find the longest carbon chain
    longest_chain_length = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:  # Carbon atom
            # Perform a depth-first search to find the longest chain starting from this atom
            visited = set()
            stack = [(atom, 1)]
            while stack:
                current_atom, current_length = stack.pop()
                visited.add(current_atom.GetIdx())
                if current_length > longest_chain_length:
                    longest_chain_length = current_length
                for neighbor in current_atom.GetNeighbors():
                    if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() not in visited:
                        stack.append((neighbor, current_length + 1))

    # Check if the longest chain length is greater than 22
    if longest_chain_length > 22:
        return True, f"Longest carbon chain length is {longest_chain_length}, which is greater than C22"
    else:
        return False, f"Longest carbon chain length is {longest_chain_length}, which is not greater than C22"