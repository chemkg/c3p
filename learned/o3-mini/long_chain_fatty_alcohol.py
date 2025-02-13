"""
Classifies: CHEBI:17135 long-chain fatty alcohol
"""
#!/usr/bin/env python
"""
Classifies: long-chain fatty alcohol
Definition: A fatty alcohol with a chain length ranging from C13 to C22.
"""

from rdkit import Chem

def is_long_chain_fatty_alcohol(smiles: str):
    """
    Determines if a molecule is a long-chain fatty alcohol.
    A long-chain fatty alcohol is defined as a compound with an aliphatic –OH group
    and a contiguous carbon chain whose length (counting only carbon atoms) falls between 13 and 22.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a long-chain fatty alcohol, False otherwise
        str: Reason for the classification
    """
    
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # First: Check for an alcohol group (–OH) attached to an aliphatic (non-aromatic) carbon.
    # We use a SMARTS query for an -OH on oxygen.
    # Note: this will also catch phenols, so we check that the neighbor carbon is not aromatic.
    alcohol_query = Chem.MolFromSmarts("[OX2H]")
    alcohol_matches = mol.GetSubstructMatches(alcohol_query)
    has_aliphatic_oh = False
    for match in alcohol_matches:
        # match is a tuple with the index of the oxygen atom
        o_idx = match[0]
        o_atom = mol.GetAtomWithIdx(o_idx)
        # Check neighbors: at least one should be a carbon that is not aromatic.
        for nbr in o_atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6 and not nbr.GetIsAromatic():
                has_aliphatic_oh = True
                break
        if has_aliphatic_oh:
            break
            
    if not has_aliphatic_oh:
        return False, "No aliphatic alcohol (-OH) functional group found"
    
    # Second: find the longest contiguous carbon chain.
    # We first gather all carbon atoms in the molecule.
    carbon_idxs = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    if not carbon_idxs:
        return False, "No carbon atoms found"
    
    # Build a dictionary mapping carbon index to a list of neighboring carbon indices.
    carbon_graph = {}
    for idx in carbon_idxs:
        atom = mol.GetAtomWithIdx(idx)
        neighbors = []
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6:  # only consider carbon neighbors
                neighbors.append(nbr.GetIdx())
        carbon_graph[idx] = neighbors
    
    # We define a DFS to find the longest simple path (without revisiting nodes) among carbon atoms.
    def dfs(node, visited):
        max_length = 1  # count the current node
        for neighbor in carbon_graph.get(node, []):
            if neighbor not in visited:
                length = 1 + dfs(neighbor, visited | {neighbor})
                if length > max_length:
                    max_length = length
        return max_length

    longest_chain = 0
    # Run DFS starting from each carbon atom.
    for idx in carbon_graph:
        chain_length = dfs(idx, {idx})
        if chain_length > longest_chain:
            longest_chain = chain_length

    # Check if the longest chain falls in the required range [13, 22].
    if 13 <= longest_chain <= 22:
        return True, f"Found an aliphatic -OH group and longest carbon chain length is {longest_chain} (within C13-C22)"
    else:
        return False, f"Longest carbon chain length is {longest_chain}, which is not in the range C13-C22"

# Example usage (uncomment for testing):
# print(is_long_chain_fatty_alcohol("CCCCCCCCCCCCCC(O)CCCC"))  # e.g., nonadecan-5-ol