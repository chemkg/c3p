"""
Classifies: CHEBI:143004 ultra-long-chain fatty acid
"""
"""
Classifies: Ultra Long‐Chain Fatty Acid (chain length > C27)
Definition: A fatty acid with a continuous, acyclic carbon chain (including the acid carbon) longer than 27 carbons.
The approach:
  1. Find free carboxylic acid group via SMARTS "[CX3](=O)[OX2H1]".
  2. From the carboxyl carbon, find the attached carbon (the fatty acid chain should be attached to the acid carbon).
  3. Perform a DFS over carbon atoms “in line” (skipping any atoms in rings) to determine the longest linear chain.
  4. If the chain length (counting the acid carbon) exceeds 27, qualify the molecule.
Note: Complex molecules might have long carbon chains in other contexts;
      we aim to restrict the search to the chain attached to the acid group.
"""

from rdkit import Chem

def is_ultra_long_chain_fatty_acid(smiles: str):
    """
    Determines whether a molecule is an ultra long-chain fatty acid,
    i.e. has a free carboxylic acid group (C(=O)[OH]) attached to a carbon chain whose
    longest continuous chain (with no ring atoms) is greater than 27.
    
    Args:
        smiles (str): Input SMILES string.
    
    Returns:
        bool: True if molecule qualifies, False otherwise.
        str: Explanation of the classification.
    """
    # Parse input SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for a free carboxylic acid group (C(=O)[OH]).
    acid_smarts = "[CX3](=O)[OX2H1]"
    acid_pattern = Chem.MolFromSmarts(acid_smarts)
    acid_matches = mol.GetSubstructMatches(acid_pattern)
    if not acid_matches:
        return False, "No free carboxylic acid group found"
    
    # Helper DFS function to find the longest acyclic carbon chain.
    # We propagate a visited set (copied at each branch) so that branching is explored.
    def dfs_chain(atom, visited):
        """
        Recursively search the longest chain (number of carbon atoms) starting from "atom".
        Only follow carbon neighbors that are not in a ring.
        """
        visited.add(atom.GetIdx())
        max_len = 1  # Count current atom
        for neighbor in atom.GetNeighbors():
            # Continue only if neighbor is carbon, has not been used, and is not in a ring.
            if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() not in visited and not neighbor.IsInRing():
                branch_len = 1 + dfs_chain(neighbor, visited.copy())
                if branch_len > max_len:
                    max_len = branch_len
        return max_len

    chain_lengths = []
    # For each acid group (if multiple exist, we consider each separately)
    for match in acid_matches:
        # In the SMARTS, index 0 is the carboxyl carbon.
        acid_carbon = mol.GetAtomWithIdx(match[0])
        # Look for carbon neighbors of the acid carbon.
        chain_neighbors = [n for n in acid_carbon.GetNeighbors() if n.GetAtomicNum() == 6]
        if not chain_neighbors:
            continue  # Not attached to any carbon chain
        # Some fatty acids might be branched at the acid carbon; choose the branch with maximal chain length.
        branch_lengths = []
        for nbr in chain_neighbors:
            # We require that the neighbor is not in a ring
            if nbr.IsInRing():
                continue
            branch_len = dfs_chain(nbr, set())
            branch_lengths.append(branch_len)
        if not branch_lengths:
            continue
        # Total chain length counts the acid carbon plus the longest branch.
        total_chain_length = 1 + max(branch_lengths)
        chain_lengths.append(total_chain_length)
    
    if not chain_lengths:
        return False, "Carboxylic acid group found but no eligible carbon chain attached"
    
    max_chain_length = max(chain_lengths)
    if max_chain_length > 27:
        return True, f"Chain length is {max_chain_length} carbons, which qualifies as an ultra-long-chain fatty acid"
    else:
        return False, f"Chain length is {max_chain_length} carbons, which does not exceed the C27 threshold"

# Example usage:
if __name__ == "__main__":
    # Test with dotriacontanoic acid: a 32-carbon fatty acid.
    test_smiles = "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC(O)=O"
    result, explanation = is_ultra_long_chain_fatty_acid(test_smiles)
    print(result, explanation)