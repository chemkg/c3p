"""
Classifies: CHEBI:17135 long-chain fatty alcohol
"""
#!/usr/bin/env python
"""
Classifies: long-chain fatty alcohol
Definition: A fatty alcohol with a chain length ranging from C13 to C22.
True long-chain fatty alcohols should have a dominant, contiguous aliphatic (non‐aromatic) carbon chain
of length between 13 and 22 carbons and at least one –OH group attached to a carbon on that chain.
False positives (such as many terpenoid, polycyclic or highly substituted structures) tend to have a long chain
but only as one fragment of a larger (often aromatic or multi‐functional) molecule.
"""

from rdkit import Chem

def is_long_chain_fatty_alcohol(smiles: str):
    """
    Determines if a molecule is a long-chain fatty alcohol.
    A long-chain fatty alcohol is defined as having an aliphatic –OH group attached 
    to a carbon that lies in a single, dominant contiguous (non-aromatic) carbon chain
    whose length is between 13 and 22 and that constitutes most of the molecule's carbons.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a long-chain fatty alcohol, False otherwise
        str: Reason for the classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # First: Check for at least one alcohol group (–OH) on an aliphatic (non‐aromatic) carbon.
    alcohol_query = Chem.MolFromSmarts("[OX2H]")  # –OH group
    alcohol_matches = mol.GetSubstructMatches(alcohol_query)
    candidate_alcohol_carbons = set()
    if not alcohol_matches:
        return False, "No -OH functional group found"
    # For each -OH oxygen found, check that it is attached to at least one carbon that is not aromatic.
    for match in alcohol_matches:
        o_idx = match[0]
        o_atom = mol.GetAtomWithIdx(o_idx)
        for nbr in o_atom.GetNeighbors():
            # if neighbor is carbon and is not aromatic then record this carbon index.
            if nbr.GetAtomicNum() == 6 and not nbr.GetIsAromatic():
                candidate_alcohol_carbons.add(nbr.GetIdx())
    if not candidate_alcohol_carbons:
        return False, "No aliphatic -OH group (attached to a non-aromatic carbon) found"
    
    # Second: Build a graph (dictionary) of connected carbon atoms (only carbons).
    carbon_idxs = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    if not carbon_idxs:
        return False, "No carbon atoms found"
        
    # Dictionary: key = carbon index, value = list of neighboring carbon indices.
    carbon_graph = {}
    for idx in carbon_idxs:
        atom = mol.GetAtomWithIdx(idx)
        neighbors = []
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6:
                neighbors.append(nbr.GetIdx())
        carbon_graph[idx] = neighbors

    # Third: Find the longest simple path (without reusing nodes) in the carbon graph.
    # We also record the actual path (list of carbon indices).
    # Because the carbon graph can have cycles, we use DFS from each node.
    longest_chain_length = 0
    longest_chain_path = []
    
    # DFS: returns the longest path (as list) starting from 'node' given a visited set.
    def dfs(node, visited):
        best_path = [node]  # current path starts with this node
        for neighbor in carbon_graph.get(node, []):
            if neighbor not in visited:
                candidate_path = [node] + dfs(neighbor, visited | {neighbor})
                if len(candidate_path) > len(best_path):
                    best_path = candidate_path
        return best_path

    for idx in carbon_graph:
        path = dfs(idx, {idx})
        if len(path) > longest_chain_length:
            longest_chain_length = len(path)
            longest_chain_path = path

    # Fourth: Check that the longest chain lies entirely on non-aromatic carbon atoms.
    for c_idx in longest_chain_path:
        if mol.GetAtomWithIdx(c_idx).GetIsAromatic():
            return False, "Longest carbon chain includes aromatic carbons"

    # Fifth: Check that the chain length is between 13 and 22.
    if not (13 <= longest_chain_length <= 22):
        return False, f"Longest carbon chain length is {longest_chain_length}, which is not in the range C13-C22"
    
    # Sixth: Check that the candidate -OH is attached to a carbon that is part of the longest chain.
    if not candidate_alcohol_carbons.intersection(longest_chain_path):
        return False, "The aliphatic -OH group is not attached to the dominant carbon chain"
    
    # Seventh: As an extra filter, require that the dominant chain represents the majority of all carbons.
    total_carbons = len(carbon_idxs)
    ratio = longest_chain_length / total_carbons
    # For “pure” fatty alcohols, nearly all carbons should be in the long chain.
    if ratio < 0.85:
        return False, (f"Longest chain (length {longest_chain_length}) accounts for only {ratio:.2f} "
                       f"of the total {total_carbons} carbons, so the molecule has significant additional fragments")
    
    return True, (f"Contains an aliphatic -OH attached to a dominant carbon chain of length "
                  f"{longest_chain_length} (which is {ratio:.2f} of all carbons) within C13-C22")

# Example usage:
# Uncomment the following lines to test a few SMILES:
# print(is_long_chain_fatty_alcohol("CCCCCCCCCCCCCCC(O)CCCC"))  # nonadecan-5-ol, expected: True
# print(is_long_chain_fatty_alcohol("C[C@@H](C1=CC=CC=C1)NC(=O)[C@@H]2[C@H](...")   # an example false positive