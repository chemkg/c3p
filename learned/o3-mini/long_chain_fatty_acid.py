"""
Classifies: CHEBI:15904 long-chain fatty acid
"""
"""
Classifies: Long-chain fatty acid (chain length ranging from C13 to C22)

A long-chain fatty acid is defined as a fatty acid with a chain length
ranging from C13 to C22. In our implementation, we require that the molecule
contains a terminal carboxylic acid group (–C(=O)O, where the acid proton might
be implicit or expressed as a negative charge) and that the longest continuous carbon
chain containing that carboxyl group (i.e. starting at its unique carbon neighbor)
has a total length (including the carboxyl carbon) between 13 and 22.
"""
from rdkit import Chem

def is_long_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a long-chain fatty acid based on its SMILES string.
    
    A long-chain fatty acid is defined as a fatty acid with a chain length ranging 
    from C13 to C22. Here we require that the molecule contains a terminal carboxyl group,
    which we (heuristically) define as the presence of a carboxyl moiety ([CX3](=O)[OX1])
    whose carbon (the acid carbon) is attached to exactly one other carbon atom. Then, 
    we compute the longest continuous carbon chain (via a DFS on the carbon-only graph)
    that is attached to that acid group. If the chain length (including the acid carbon)
    is between 13 and 22 (inclusive), we return True.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule qualifies as a long-chain fatty acid, False otherwise.
        str: Reason for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # We define a SMARTS for a carboxyl group.
    # This pattern [CX3](=O)[OX1] will match a carbonyl carbon attached to an oxygen,
    # either as -OH or as [O-]. This is more general than requiring an explicit H.
    acid_smarts = "[CX3](=O)[OX1]"
    acid_pattern = Chem.MolFromSmarts(acid_smarts)
    acid_matches = mol.GetSubstructMatches(acid_pattern)
    if not acid_matches:
        return False, "Molecule does not have a carboxyl group (C(=O)O or C(=O)[O-])"
    
    # Build a carbon connectivity graph for the molecule.
    # The graph keys are atom indices of carbon atoms (atomic number 6)
    # and the values are lists of bonded carbon atom indices.
    carbon_graph = {}
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:
            idx = atom.GetIdx()
            carbon_graph[idx] = []
            for nbr in atom.GetNeighbors():
                if nbr.GetAtomicNum() == 6:
                    carbon_graph[idx].append(nbr.GetIdx())
    
    # Now search among all acid matches for one that is "terminal".
    # We require that the acid carbon (first atom in the match) has exactly one bonded carbon neighbor.
    valid_acid_found = False
    chain_length_found = None
    acid_atom_idx = None
    chain_start_idx = None
    for match in acid_matches:
        # In the SMARTS pattern "[CX3](=O)[OX1]", the first atom is the acid carbon.
        acid_idx = match[0]
        # Get the carbon neighbors of this acid carbon
        acid_atom = mol.GetAtomWithIdx(acid_idx)
        carbon_neighbors = [nbr.GetIdx() for nbr in acid_atom.GetNeighbors() if nbr.GetAtomicNum() == 6]
        if len(carbon_neighbors) != 1:
            # Not terminal (could be attached to more than one carbon); skip this match.
            continue
        # We have found a candidate. Record the acid carbon and its unique carbon neighbor.
        acid_atom_idx = acid_idx
        chain_start_idx = carbon_neighbors[0]
        valid_acid_found = True
        break

    if not valid_acid_found:
        return False, "No terminal carboxyl group found (the acid carbon is not at the chain terminus)"
    
    # Now compute the length of the long carbon chain that is connected to the acid group.
    # We define a recursive DFS to find the longest simple path in the carbon_graph.
    def dfs(current, visited):
        max_length = 1  # Count the current node
        for nbr in carbon_graph.get(current, []):
            if nbr in visited:
                continue
            new_visited = visited | {nbr}
            candidate = 1 + dfs(nbr, new_visited)
            if candidate > max_length:
                max_length = candidate
        return max_length

    # Our chain includes the acid carbon plus the chain starting from its unique carbon neighbor.
    # Start the DFS from the chain_start_idx; mark the acid carbon as visited so we don't go back.
    chain_length = 1 + dfs(chain_start_idx, {acid_atom_idx, chain_start_idx})
    
    # Check if the measured chain length falls within the allowed range
    if chain_length < 13:
        return False, f"Chain length too short: {chain_length} carbons (< 13 required)"
    if chain_length > 22:
        return False, f"Chain length too long: {chain_length} carbons (> 22 allowed)"
    
    return True, f"Terminal carboxyl group found with chain length of {chain_length} carbons"

# Example test cases (uncomment for local testing):
# examples = [
#     ("O(O)[C@H](CCCCC)\\C=C\\CCCCCCCCCC(O)=O", "13R-HpOME(11E)"),
#     ("O[C@H]([C@@H](O)C/C=C/CCCCC)C(=O)CCCCCCC(O)=O", "(9R,10S,12Z)-9,10-Dihydroxy-8-oxo-12-octadecenoic acid"),
#     ("CCCCCCCCCCCCC(O)=O", "Tridecanoic acid"),
# ]
# for s, name in examples:
#     result, reason = is_long_chain_fatty_acid(s)
#     print(f"SMILES: {s}\nNAME: {name}\nResult: {result}, Reason: {reason}\n")
    
# Note: This heuristic defines the chain based on the longest simple path through carbon atoms
# that is connected at the terminal carboxyl group. This approach may need refinement for
# very complex structures, but it addresses the issues encountered previously.