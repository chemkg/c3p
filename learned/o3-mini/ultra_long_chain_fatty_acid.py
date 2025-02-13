"""
Classifies: CHEBI:143004 ultra-long-chain fatty acid
"""
"""
Classifies: Ultra Long-Chain Fatty Acid (chain length > C27)
Definition: A fatty acid which has a continuous carbon chain (including the acid carbon) longer than 27 carbons.
"""
from rdkit import Chem

def is_ultra_long_chain_fatty_acid(smiles: str):
    """
    Determines whether a molecule is an ultra long-chain fatty acid based on its SMILES.
    It checks for a carboxylic acid group and then finds the longest continuous carbon chain 
    connected to that group. The chain length (number of carbons, including the carboxyl carbon)
    must be greater than 27.
    
    Args:
        smiles (str): Input SMILES string for the molecule.
        
    Returns:
        bool: True if the molecule qualifies as an ultra-long-chain fatty acid, False otherwise.
        str: Explanation of the classification.
    """
    # Parse the input SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a simple SMARTS pattern to identify a carboxylic acid group.
    # (Note: This might catch ester groups as well, but in typical fatty acids the acid group is free.)
    acid_pattern = Chem.MolFromSmarts("C(=O)O")
    acid_matches = mol.GetSubstructMatches(acid_pattern)
    if not acid_matches:
        return False, "No carboxylic acid group found"
    
    # We define a helper function to perform a depth-first search (DFS) along carbon atoms.
    def dfs(atom, parent, visited):
        """
        Recursively computes the length of the carbon chain starting from 'atom'.
        We count 'atom' itself and then traverse to adjacent carbon atoms (atomic number 6)
        that have not yet been visited.
        """
        visited.add(atom.GetIdx())
        max_length = 1  # Count the current atom.
        for neighbor in atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() != parent.GetIdx() and neighbor.GetIdx() not in visited:
                branch_length = 1 + dfs(neighbor, atom, visited)
                if branch_length > max_length:
                    max_length = branch_length
        return max_length

    chain_lengths = []
    # Loop over each matched carboxylic acid group.
    for match in acid_matches:
        # The first atom in the matched SMARTS (index 0 in the match tuple) corresponds to the carboxyl carbon.
        acid_carbon = mol.GetAtomWithIdx(match[0])
        # In a typical fatty acid the acid carbon will be connected to two oxygens and exactly one carbon.
        # So we look for that one carbon neighbor.
        chain_neighbors = [n for n in acid_carbon.GetNeighbors() if n.GetAtomicNum() == 6]
        if not chain_neighbors:
            continue  # Skip if no carbon neighbor exists.
        # Use a new visited set for each acid group search.
        visited = set()
        # We start the DFS from the acid carbon. The chain length here includes the acid carbon itself.
        total_chain_length = dfs(acid_carbon, acid_carbon, visited)
        chain_lengths.append(total_chain_length)
    
    if not chain_lengths:
        return False, "Carboxylic acid group found but no attached carbon chain detected"
    
    # Use the maximum chain length from any identified acid group.
    max_chain_length = max(chain_lengths)
    # Check if the chain length (number of carbons) exceeds 27.
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