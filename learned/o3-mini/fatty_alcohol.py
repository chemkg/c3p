"""
Classifies: CHEBI:24026 fatty alcohol
"""
"""
Classifies: Fatty alcohol
Definition: An aliphatic alcohol consisting of a chain of 3 to greater than 27 carbon atoms.
Fatty alcohols may be saturated or unsaturated and may be branched or unbranched.
This implementation adds explicit hydrogens and then identifies candidate -OH groups
on non‐aromatic carbons not involved in carbonyls. For each candidate, it builds a graph
of connected non‐aromatic carbons to compute the longest contiguous chain attached.
The molecule is accepted if at least one candidate –OH connects to a contiguous chain of at least 3 carbons.
"""

from rdkit import Chem
from collections import deque

def is_fatty_alcohol(smiles: str):
    """
    Determines if a molecule is a fatty alcohol based on its SMILES string.
    A fatty alcohol is defined as an aliphatic alcohol having a contiguous chain (of non‐aromatic
    carbons) of at least 3 atoms attached to a valid hydroxyl (-OH) group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty alcohol, False otherwise
        str: Reason for classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    # Add explicit hydrogens so that hydroxyl groups can be correctly detected.
    mol = Chem.AddHs(mol)
    
    # Step 1: Detect candidate aliphatic hydroxyl groups.
    # A valid candidate oxygen should have exactly one hydrogen neighbor and at least one non‐aromatic carbon neighbor.
    # Also, the carbon neighbor should not be in a carbonyl context (i.e. double-bonded to another oxygen).
    candidate_carbon_idxs = set()
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 8:  # oxygen
            # Count hydrogen neighbors
            h_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 1]
            if len(h_neighbors) != 1:
                continue
            # Look for a connected carbon that is non-aromatic and not in a carbonyl environment.
            for nbr in atom.GetNeighbors():
                if nbr.GetAtomicNum() == 6 and not nbr.GetIsAromatic():
                    is_carbonyl = False
                    # Check bonds from this carbon: if any oxygen neighbor bond is a double bond, then skip.
                    for cnbr in nbr.GetNeighbors():
                        if cnbr.GetAtomicNum() == 8 and cnbr.GetIdx() != atom.GetIdx():
                            bond = mol.GetBondBetweenAtoms(nbr.GetIdx(), cnbr.GetIdx())
                            if bond is not None and bond.GetBondType() == Chem.BondType.DOUBLE:
                                is_carbonyl = True
                                break
                    if not is_carbonyl:
                        candidate_carbon_idxs.add(nbr.GetIdx())
    
    if not candidate_carbon_idxs:
        return False, "No valid aliphatic hydroxyl group found (or all -OH groups are in non-alcoholic contexts)"
    
    # Step 2: Build a graph of all non-aromatic carbon atoms.
    # The graph is represented as a dictionary that maps atom indices to a set of connected carbon atom indices.
    carbon_graph = {}
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6 and not atom.GetIsAromatic():
            carbon_graph[atom.GetIdx()] = set()
    
    # Add edges for bonds between these non-aromatic carbons.
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        if (a1.GetAtomicNum() == 6 and not a1.GetIsAromatic() and
            a2.GetAtomicNum() == 6 and not a2.GetIsAromatic()):
            idx1 = a1.GetIdx()
            idx2 = a2.GetIdx()
            if idx1 in carbon_graph and idx2 in carbon_graph:
                carbon_graph[idx1].add(idx2)
                carbon_graph[idx2].add(idx1)
    
    if not carbon_graph:
        return False, "No non-aromatic carbon chain found"
    
    # Helper function: perform a BFS on an undirected graph starting from a given node.
    # Returns the farthest node and a dictionary mapping each node to its distance from the start.
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
        # Determine the farthest node and its distance.
        farthest_node = start
        max_distance = 0
        for node, dist in distances.items():
            if dist > max_distance:
                max_distance = dist
                farthest_node = node
        return farthest_node, distances

    # Helper function: given a connected component (set of nodes) in the carbon graph,
    # compute its diameter (the longest shortest-path in the component). The chain length in atoms is diameter + 1.
    def component_diameter(component_nodes, graph):
        start = next(iter(component_nodes))
        node_a, _ = bfs(start, graph)
        node_b, distances = bfs(node_a, graph)
        return max(distances.values()) + 1

    # Step 3: For each candidate carbon (the one attached to a valid hydroxyl),
    # find the connected component (in the non-aromatic carbon graph)
    # and compute the longest chain (its diameter) in that component.
    best_chain = 0
    processed_components = set()  # to avoid repeating work on the same component
    
    for candidate in candidate_carbon_idxs:
        if candidate not in carbon_graph:
            continue
        # If we already processed this component, skip it.
        if candidate in processed_components:
            continue
        
        # Find the full connected component by DFS.
        component_nodes = set()
        stack = [candidate]
        while stack:
            node = stack.pop()
            if node in component_nodes:
                continue
            component_nodes.add(node)
            for nbr in carbon_graph.get(node, []):
                if nbr not in component_nodes:
                    stack.append(nbr)
        
        processed_components.update(component_nodes)
        chain_length = component_diameter(component_nodes, carbon_graph)
        # We require a chain of at least 3 carbon atoms.
        if chain_length >= 3 and chain_length > best_chain:
            best_chain = chain_length
    
    if best_chain >= 3:
        return True, f"Found a contiguous aliphatic carbon chain of {best_chain} atoms attached to a valid hydroxyl group."
    else:
        return False, "No contiguous aliphatic carbon chain of minimum length (3) found attached to a valid hydroxyl group."

# For testing (this section can be removed in production code)
if __name__ == "__main__":
    test_smiles = [
        "CCCCCCCCCCCCCCC(O)CCCC",  # nonadecan-5-ol
        "O(C(CCCCCCCCCC(O)C)CC1=CC(O)=CC(O)=C1)C(=O)C",  # [1-(3,5-Dihydroxyphenyl)-12-hydroxytridecan-2-yl] acetate
        "CCCCCCCCCCCCCC(O)CCCC",  # octadecan-5-ol
        "OC(CCCCCC)CC(=O)CCCCCCCCCCCCC",  # 6-Hydroxy-8-heneicosanone
        "CCCCCCCCCCCCCCC(O)CCCCCCCCCC",  # pentacosan-11-ol
        "OC(C(O)CO)C(O)C#CC#CC#CC#CC",  # (2S,3S,4S)-5,7,9,11-Tridecatetrayne-1,2,3,4-tetrol
        "O[C@H](/C(=C/[C@H]([C@H](O)C)C)/C)C/C=C/C=C/C=C/C[C@@H](O)[C@@H](O)CC",  # Succinilene D
        "CCCCCCCC(O)CCC",  # undecan-4-ol
        "CCCCCCCCCCC(O)CCCCC",  # hexadecan-6-ol
        "CCCCCCC(O)CCCCCCCO",  # 1,8-tetradecanediol
    ]
    for s in test_smiles:
        result, reason = is_fatty_alcohol(s)
        print(f"SMILES: {s}\nResult: {result}\nReason: {reason}\n{'-'*60}")