"""
Classifies: CHEBI:73702 wax
"""
"""
Classifies: Wax
Definition: A wax is an organic compound or mixture of compounds composed of long‐chain molecules and is malleable at ambient temperatures.
This classifier uses the following heuristics:
  - The molecule must parse correctly.
  - The molecular weight must be above 300 Da.
  - The molecule must have at least 15 carbon atoms.
  - The molecule must be mostly aliphatic — if more than 30% of its carbon atoms are aromatic, we reject it.
  - The molecule must contain at least one uninterrupted carbon chain of 12 or more carbons.
  
These heuristics help avoid costly computations on highly cyclic (eg, aromatic) systems that are unlikely to be waxes.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from functools import lru_cache

def longest_carbon_chain_length(mol):
    """
    Calculate the longest continuous chain of carbon atoms (atomic number 6).
    Uses a DFS on the subgraph containing only carbon atoms.
    Memoization is applied to avoid redundant calculations.
    
    Returns:
        int: The number of carbon atoms in the longest chain.
    """
    # Get indices of carbon atoms
    carbon_indices = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    if not carbon_indices:
        return 0

    # Build an adjacency list for carbon atoms only
    carbon_neighbors = {idx: [] for idx in carbon_indices}
    for idx in carbon_indices:
        atom = mol.GetAtomWithIdx(idx)
        for neighbor in atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 6:
                carbon_neighbors[idx].append(neighbor.GetIdx())

    # Use memoization in DFS to avoid repeated work.
    @lru_cache(maxsize=None)
    def dfs(current, visited):
        max_length = 1  # current atom counts as 1
        for neighbor in carbon_neighbors[current]:
            if neighbor not in visited:
                new_visited = visited + (neighbor,)  # using tuple to be hashable
                new_length = 1 + dfs(neighbor, new_visited)
                if new_length > max_length:
                    max_length = new_length
        return max_length

    longest = 0
    # Try every carbon atom as starting point
    for start in carbon_indices:
        chain_length = dfs(start, (start,))
        if chain_length > longest:
            longest = chain_length
    return longest

def is_wax(smiles: str):
    """
    Classifies a molecule as a wax based on its SMILES string.
    
    Heuristics used:
      - Successful parsing.
      - Molecular weight > 300 Da.
      - At least 15 carbon atoms.
      - Majority aliphatic: less than 30% aromatic carbon atoms.
      - Contains at least one uninterrupted carbon chain of 12 carbons or more.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if classified as a wax, False otherwise.
        str: Explanation for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # Check molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:
        return False, f"Molecular weight too low for a wax ({mol_wt:.1f} Da)."
    
    # Count carbon atoms and aromatic carbon atoms
    carbon_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    num_carbons = len(carbon_atoms)
    if num_carbons < 15:
        return False, f"Insufficient carbon count ({num_carbons}). Must have at least 15 carbons for a wax."
    
    aromatic_carbons = [atom for atom in carbon_atoms if atom.GetIsAromatic()]
    frac_aromatic = len(aromatic_carbons) / num_carbons
    if frac_aromatic >= 0.3:
        return False, f"Too many aromatic carbons ({len(aromatic_carbons)} out of {num_carbons}) - likely not a wax."
    
    # Determine the length of the longest uninterrupted carbon chain
    longest_chain = longest_carbon_chain_length(mol)
    if longest_chain < 12:
        return False, f"Longest continuous carbon chain is only {longest_chain} atoms (requires at least 12)."
    
    return True, f"Valid wax: Molecular weight {mol_wt:.1f} Da, {num_carbons} carbons (aromatic fraction {frac_aromatic:.2f}), longest chain of {longest_chain} carbons."

# Example usage (uncomment the lines below to test):
# wax_examples = [
#     "O(CCCCCCCC/C=C\\C/C=C\\CCCCC)C(=O)CCCCCCC/C=C\\C/C=C\\CCCCC",  # Linoleyl linoleate
#     "CCCCCCCCCCCCCCCCCCOC(=O)CCCCCCC/C=C\\CCCCCCCC",  # stearyl oleate (variation example)
# ]
# for s in wax_examples:
#     result, reason = is_wax(s)
#     print(s)
#     print(result, reason)