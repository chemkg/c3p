"""
Classifies: CHEBI:24026 fatty alcohol
"""
"""
Classifies: Fatty alcohol
Definition: An aliphatic alcohol consisting of a chain of 3 to greater than 27 carbon atoms. 
Fatty alcohols may be saturated or unsaturated and may be branched or unbranched.
Improvements over the previous attempt:
  - Only consider hydroxyl groups (-OH) on non‐aromatic carbons that are not part of carboxyl (or carbonyl) groups.
  - Build a graph of all connected non‐aromatic carbon atoms and compute the longest chain (diameter) in that graph.
  - Accept the molecule if at least one candidate –OH connects to a contiguous chain of at least 3 carbon atoms.
"""
from rdkit import Chem
from collections import deque

def is_fatty_alcohol(smiles: str):
    """
    Determines if a molecule is a fatty alcohol based on its SMILES string.
    A fatty alcohol is defined as an aliphatic alcohol having a chain (of non‐aromatic
    carbons) of at least 3 atoms attached to a valid hydroxyl group. (The chain can be saturated,
    unsaturated, and may be branched.)
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a fatty alcohol, False otherwise.
        str: Reason for classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Step 1: Identify candidate (aliphatic) hydroxyl groups.
    # A valid candidate is an oxygen (atomic num 8) having exactly one hydrogen neighbor,
    # attached to at least one carbon atom that is non-aromatic and is not part of an acid/carbonyl group.
    candidate_idxs = set()
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 8:  # oxygen
            # Check that it has exactly one hydrogen neighbor (an -OH group)
            h_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 1]
            if len(h_neighbors) != 1:
                continue
            # Look for a carbon neighbor
            for nbr in atom.GetNeighbors():
                if nbr.GetAtomicNum() == 6 and not nbr.GetIsAromatic():
                    # Exclude if this carbon is in a carboxyl or carbonyl group:
                    # i.e. if any oxygen neighbor is double bonded to this carbon.
                    is_carboxyl = False
                    for cnbr in nbr.GetNeighbors():
                        if cnbr.GetAtomicNum() == 8:
                            # Look at the bond between nbr and cnbr
                            bond = mol.GetBondBetweenAtoms(nbr.GetIdx(), cnbr.GetIdx())
                            if bond is not None and bond.GetBondType() == Chem.BondType.DOUBLE:
                                is_carboxyl = True
                                break
                    if not is_carboxyl:
                        candidate_idxs.add(nbr.GetIdx())
    if not candidate_idxs:
        return False, "No valid aliphatic hydroxyl group found (or all -OH groups are in non-alcoholic contexts)"
    
    # Step 2: Build a graph (as a dictionary) of all non-aromatic carbon atoms.
    # Each node key is the atom index and value is a set of neighbor indices (only carbons).
    carbon_graph = {}
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6 and not atom.GetIsAromatic():
            idx = atom.GetIdx()
            carbon_graph[idx] = set()
    # Now add edges (for bonds between these carbons)
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        if (a1.GetAtomicNum() == 6 and (not a1.GetIsAromatic()) and 
            a2.GetAtomicNum() == 6 and (not a2.GetIsAromatic())):
            idx1 = a1.GetIdx()
            idx2 = a2.GetIdx()
            if idx1 in carbon_graph and idx2 in carbon_graph:
                carbon_graph[idx1].add(idx2)
                carbon_graph[idx2].add(idx1)
    
    if not carbon_graph:
        return False, "No non-aromatic carbon chain found"
    
    # Helper function: perform BFS in an undirected graph from a starting node.
    # Returns the farthest node and a dictionary mapping node->distance.
    def bfs(start, graph):
        visited = {start}
        distances = {start: 0}
        queue = deque([start])
        while queue:
            current = queue.popleft()
            for neighbor in graph.get(current, []):
                if neighbor not in visited:
                    visited.add(neighbor)
                    distances[neighbor] = distances[current] + 1
                    queue.append(neighbor)
        # Find the farthest node and its distance.
        farthest_node = start
        max_distance = 0
        for node, dist in distances.items():
            if dist > max_distance:
                max_distance = dist
                farthest_node = node
        return farthest_node, distances

    # Helper function: given a set of nodes that form a connected component (subgraph),
    # compute its diameter (longest shortest-path between any two nodes) using two BFS runs.
    def component_diameter(component_nodes, graph):
        # Pick an arbitrary node from the component.
        start = next(iter(component_nodes))
        node_a, _ = bfs(start, graph)
        node_b, distances = bfs(node_a, graph)
        # The longest distance (in number of bonds) found is the diameter (number of edges)
        # The chain length in atoms is diameter + 1.
        return max(distances.values()) + 1

    # Step 3: For each candidate carbon (attached to a valid hydroxyl), compute the connected component
    # (in the non-aromatic carbon graph) that contains it, and then compute the longest chain (diameter) of that component.
    best_chain = 0
    best_candidate = None
    visited_components = set()  # to avoid recomputing for the same component
    for candidate in candidate_idxs:
        if candidate not in carbon_graph:
            # In some cases the candidate might not be in our carbon graph (if isolated etc.)
            continue
        if candidate in visited_components:
            # Already computed for this component.
            continue
        # Get all nodes in the connected component
        comp_nodes = set()
        stack = [candidate]
        while stack:
            node = stack.pop()
            if node in comp_nodes:
                continue
            comp_nodes.add(node)
            for nbr in carbon_graph.get(node, []):
                if nbr not in comp_nodes:
                    stack.append(nbr)
        # Mark these nodes as visited so we do not reprocess them
        visited_components.update(comp_nodes)
        # Compute the diameter of this component as the longest aliphatic chain in it
        chain_length = component_diameter(comp_nodes, carbon_graph)
        if chain_length > best_chain:
            best_chain = chain_length
            best_candidate = candidate

    # Step 4: Decide if the molecule is a fatty alcohol based on the longest chain found.
    # According to the definition the chain should be at least 3 carbons long.
    if best_chain >= 3:
        return True, f"Found a contiguous aliphatic carbon chain of {best_chain} atoms attached to a valid hydroxyl group."
    else:
        return False, "No contiguous aliphatic carbon chain of minimum length (3) found attached to a valid hydroxyl group."

# For testing purposes (can be removed in a production module):
if __name__ == "__main__":
    test_smiles = [
        "CCCCCCCCCCCCCCC(O)CCCC",  # nonadecan-5-ol: should be classified as fatty alcohol.
        "C(C(C([2H])([2H])[2H])(N([2H])[2H])[2H])(=O)O[2H]",  # alanine-d7: should NOT be classified.
        "CCCCCCCCCCCCCC(O)CCCC",  # octadecan-5-ol, etc.
    ]
    for s in test_smiles:
        result, reason = is_fatty_alcohol(s)
        print(f"SMILES: {s}\nResult: {result}\nReason: {reason}\n{'-'*60}")