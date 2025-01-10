"""
Classifies: CHEBI:143004 ultra-long-chain fatty acid
"""
from rdkit import Chem

def is_ultra_long_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is an ultra-long-chain fatty acid based on its SMILES string.
    Defined as a fatty acid with a chain length greater than 27 carbon atoms.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an ultra-long-chain fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Identify carboxylic acid group or similar termini
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)O")
    matches = mol.GetSubstructMatches(carboxyl_pattern)
    if not matches:
        return False, "No carboxyl group found"
        
    max_chain_length = 0
    
    # Use Breadth-First Search to measure chain lengths from each carbon
    for match in matches:
        carboxyl_carbon = match[0]
        visited = set()
        queue = [(carboxyl_carbon, 0)]
        
        while queue:
            atom_idx, chain_length = queue.pop(0)
            if atom_idx in visited:
                continue
            visited.add(atom_idx)
            
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetSymbol() == 'C':
                chain_length += 1
                for neighbor in atom.GetNeighbors():
                    neighbor_idx = neighbor.GetIdx()
                    if neighbor_idx not in visited and neighbor.GetSymbol() == 'C':
                        queue.append((neighbor_idx, chain_length))
            
            max_chain_length = max(max_chain_length, chain_length)
    
    if max_chain_length > 27:
        return True, f"Contains {max_chain_length} carbon atoms in chain, qualifies as ultra-long-chain fatty acid"
    else:
        return False, f"Contains {max_chain_length} carbon atoms in chain, does not qualify as ultra-long-chain fatty acid"