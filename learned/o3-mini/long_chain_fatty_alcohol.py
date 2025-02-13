"""
Classifies: CHEBI:17135 long-chain fatty alcohol
"""
#!/usr/bin/env python
"""
Classifies: long-chain fatty alcohol
Definition: A fatty alcohol with a chain length ranging from C13 to C22.
True long-chain fatty alcohols should have an aliphatic –OH group attached to a carbon that is part of a
single, dominant, contiguous (non‐ring, non‐aromatic) aliphatic carbon chain whose length lies between 13 and 22
and which constitutes most of the molecule's carbons.
"""

from rdkit import Chem

def is_long_chain_fatty_alcohol(smiles: str):
    """
    Determines if a molecule is a long-chain fatty alcohol.
    A long-chain fatty alcohol is defined as having an aliphatic –OH group attached to a carbon (which is in
    a dominant, contiguous aliphatic chain) and the longest such chain (composed solely of carbons that are
    not aromatic and not in any ring) has a length between 13 and 22. In addition, the dominant chain should
    represent most of the carbons in the molecule.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a long-chain fatty alcohol, False otherwise.
        str: Reason for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Step 1: Identify candidate -OH groups.
    # We look for an –OH and record those carbons that are non-aromatic neighbors.
    alcohol_query = Chem.MolFromSmarts("[OX2H]")  # –OH group
    alcohol_matches = mol.GetSubstructMatches(alcohol_query)
    candidate_alcohol_carbons = set()
    if not alcohol_matches:
        return False, "No -OH functional group found"
    
    for match in alcohol_matches:
        o_idx = match[0]
        o_atom = mol.GetAtomWithIdx(o_idx)
        for nbr in o_atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6 and not nbr.GetIsAromatic():
                candidate_alcohol_carbons.add(nbr.GetIdx())
    
    if not candidate_alcohol_carbons:
        return False, "No aliphatic -OH group (attached to a non-aromatic carbon) found"
    
    # Step 2: Determine the total carbon count (all carbons in the molecule).
    total_carbons = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    if not total_carbons:
        return False, "No carbon atoms found"
    
    # Step 3: Identify all carbon atoms that are in rings.
    ring_info = mol.GetRingInfo()
    ring_carbons = set()
    for ring in ring_info.AtomRings():
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            # We add the carbon regardless whether it is aromatic or not;
            # for our definition of a fatty alcohol we want a contiguous chain from non‐ring atoms.
            if atom.GetAtomicNum() == 6:
                ring_carbons.add(idx)
    
    # Step 4: Create a set of candidate carbon atoms in the acyclic aliphatic part.
    # Conditions: must be carbon, not aromatic and not part of any ring.
    nonring_aliphatic = set()
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6 and (not atom.GetIsAromatic()) and atom.GetIdx() not in ring_carbons:
            nonring_aliphatic.add(atom.GetIdx())
    
    if not nonring_aliphatic:
        return False, "No acyclic aliphatic carbon chain found (all carbons are in rings or aromatic)"
    
    # Step 5: Build a graph (dictionary) of connected acyclic aliphatic carbons.
    carbon_graph = {}
    for idx in nonring_aliphatic:
        atom = mol.GetAtomWithIdx(idx)
        neighbors = []
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6 and nbr.GetIdx() in nonring_aliphatic:
                neighbors.append(nbr.GetIdx())
        carbon_graph[idx] = neighbors

    # Step 6: Find the longest simple path in the graph.
    longest_chain_length = 0
    longest_chain_path = []

    def dfs(node, visited):
        # DFS starting at node, return the longest simple path (list of nodes).
        best_path = [node]
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

    if longest_chain_length == 0:
        return False, "Could not determine a contiguous aliphatic carbon chain"
    
    # Step 7: Ensure that the longest chain (by non-ring definition) is within C13-C22.
    if not (13 <= longest_chain_length <= 22):
        return False, (f"Longest carbon chain length is {longest_chain_length}, "
                       f"which is not in the required range C13-C22")
    
    # Step 8: Ensure the candidate -OH is attached to a carbon that is part of the longest chain.
    if not candidate_alcohol_carbons.intersection(longest_chain_path):
        return False, "The aliphatic -OH group is not attached to the dominant carbon chain"
    
    # Step 9: Check that the dominant chain accounts for the majority of the molecule's carbon atoms.
    total_carbon_count = len(total_carbons)
    ratio = longest_chain_length / total_carbon_count
    # We require that at least 85% of the molecule’s carbons (or a very high fraction) belong to the chain.
    if ratio < 0.85:
        return False, (f"Longest chain (length {longest_chain_length}) accounts for only {ratio:.2f} "
                       f"of the total {total_carbon_count} carbons, so the molecule has significant additional fragments")
    
    return True, (f"Contains an aliphatic -OH attached to a dominant carbon chain of length "
                  f"{longest_chain_length} (which is {ratio:.2f} of all carbons) within C13-C22")

# Example usage:
# To test a few SMILES strings, you can uncomment the following lines:
# print(is_long_chain_fatty_alcohol("CCCCCCCCCCCCCC(O)CCCC"))  # nonadecan-5-ol, expected True
# print(is_long_chain_fatty_alcohol("O1[C@@H]2[C@H](O)C(/C=C/[C@H](O)CCCCCCCC)=C([C@H]([C@H]12)O)CO")  # Phomopoxide D, expected True