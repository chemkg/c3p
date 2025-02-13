"""
Classifies: CHEBI:73702 wax
"""
"""
Classifies: Wax
Definition: A chemical substance that is an organic compound or mixture of compounds that is composed of long-chain molecules and is malleable at ambient temperatures.
This simple classifier uses a set of heuristics:
    - The molecule must be valid.
    - The molecule must have a high enough molecular weight (>300 Da) and a sufficient number of carbons.
    - The molecule must feature at least one long uninterrupted carbon chain (>= 12 carbons in a row).
Many waxes (for example, fatty acid esters) fulfill these criteria.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def longest_carbon_chain_length(mol):
    """
    Calculate the longest continuous chain of carbon atoms (atomic number 6).
    Uses a depth-first search (DFS) on the subgraph consisting only of carbons.
    
    Returns:
        int: The number of carbon atoms in the longest chain.
    """
    # Build a mapping: atom index -> atom object if carbon
    carbon_indices = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    if not carbon_indices:
        return 0

    # Create an adjacency list for carbon atoms only
    carbon_neighbors = {idx: [] for idx in carbon_indices}
    for idx in carbon_indices:
        atom = mol.GetAtomWithIdx(idx)
        for neighbor in atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 6:
                carbon_neighbors[idx].append(neighbor.GetIdx())

    # Depth-first search (DFS) routine to find the longest simple path
    def dfs(current, visited):
        max_length = 1  # current atom counts as 1
        for neighbor in carbon_neighbors[current]:
            if neighbor not in visited:
                # Visit neighbor and update path length
                new_length = 1 + dfs(neighbor, visited | {neighbor})
                if new_length > max_length:
                    max_length = new_length
        return max_length

    longest_chain = 0
    for start in carbon_indices:
        chain_len = dfs(start, {start})
        if chain_len > longest_chain:
            longest_chain = chain_len
    return longest_chain

def is_wax(smiles: str):
    """
    Classifies a molecule as a wax based on its SMILES string.
    
    Heuristics used:
      - The molecule must be parsed successfully.
      - Its molecular weight should be high enough (>300 Da).
      - It must contain a sufficient number of carbon atoms (>=15) to be considered "long-chain".
      - It must have at least one uninterrupted carbon chain of at least 12 carbon atoms.
      
    These heuristics were chosen based on the definition of waxes and examples of long-chain fatty acid esters.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if the molecule is classified as a wax, False otherwise.
        str: A reason explaining the result.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # Calculate molecular weight.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:
        return False, f"Molecular weight too low for a wax ({mol_wt:.1f} Da)."
    
    # Count the total number of carbon atoms.
    c_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    if len(c_atoms) < 15:
        return False, f"Not enough carbon atoms ({len(c_atoms)}), expected at least 15 for a wax."

    # Determine the longest uninterrupted carbon chain.
    longest_chain = longest_carbon_chain_length(mol)
    if longest_chain < 12:
        return False, f"Longest continuous carbon chain is only {longest_chain} atoms (requires at least 12)."
        
    # If all heuristics pass, we classify molecule as a wax.
    return True, f"Valid wax: Molecular weight {mol_wt:.1f} Da, {len(c_atoms)} carbons, longest chain {longest_chain}."

# Example usage (uncomment to test):
# smiles_examples = [
#     "O(CCCCCCCC/C=C\\C/C=C\\CCCCC)C(=O)CCCCCCC/C=C\\C/C=C\\CCCCC",  # Linoleyl linoleate
#     "CCCCCCCCCCCCCCCCCCOC(=O)CCCCCCC/C=C\\CCCCCCCC",  # stearyl oleate (example variation)
# ]
# for s in smiles_examples:
#     result, reason = is_wax(s)
#     print(s, result, reason)