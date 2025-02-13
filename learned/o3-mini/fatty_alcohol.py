"""
Classifies: CHEBI:24026 fatty alcohol
"""
"""
Classifies: Fatty alcohol
Definition: An aliphatic alcohol consisting of a chain of 3 (for small molecules) to >27 carbon atoms.
Fatty alcohols may be saturated or unsaturated and may be branched or unbranched.
This implementation:
  1. Adds explicit hydrogens.
  2. Identifies candidate –OH groups (on non‐aromatic carbons that are not in a carbonyl).
  3. Builds a graph from all non–aromatic carbon atoms.
  4. For candidates, finds the connected carbon component and computes its “diameter” (chain length).
  5. Applies two extra “filters”:
       a. For larger molecules (with many heavy atoms) the chain length must be at least 7:
          (small molecules, e.g. propanol, are allowed if total heavy atoms ≤20 and chain length ≥3).
       b. If more than two –OH groups (in “simple” aliphatic contexts) are present the molecule is rejected.
The extra rules are ad hoc choices aimed at reducing false positives.
"""

from rdkit import Chem
from collections import deque

def is_fatty_alcohol(smiles: str):
    """
    Determines if a molecule is a fatty alcohol based on its SMILES string.
    A fatty alcohol is an aliphatic alcohol in which at least one valid –OH (on a non‐aromatic,
    non–carbonyl carbon) is attached to a contiguous chain of non‐aromatic carbons.
    For small molecules (≤20 heavy atoms) a chain of at least 3 carbons qualifies,
    while for larger molecules the chain must be at least 7 carbons long.
    In addition, if the molecule has more than 2 such hydroxyl groups it is presumed too complex.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is deemed a fatty alcohol, False otherwise.
        str: Reason for classification.
    """
    # Parse SMILES and add explicit hydrogens.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    mol = Chem.AddHs(mol)
    
    # Step 1: Identify candidate –OH groups.
    # We require an oxygen atom with exactly one H neighbor and attached to at least one
    # non‐aromatic carbon that is not in a carbonyl environment.
    candidate_carbon_idxs = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 8:
            continue
        # Count hydrogen neighbors (explicit H)
        h_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 1]
        if len(h_neighbors) != 1:
            continue
        # Look for a connected carbon (non–aromatic) that is not in a carbonyl.
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() != 6 or nbr.GetIsAromatic():
                continue
            is_carbonyl = False
            for cnbr in nbr.GetNeighbors():
                # if a carbon is bonded (via a double bond) to another oxygen (other than the -OH)
                if cnbr.GetAtomicNum() == 8 and cnbr.GetIdx() != atom.GetIdx():
                    bond = mol.GetBondBetweenAtoms(nbr.GetIdx(), cnbr.GetIdx())
                    if bond is not None and bond.GetBondType() == Chem.BondType.DOUBLE:
                        is_carbonyl = True
                        break
            if not is_carbonyl:
                candidate_carbon_idxs.append(nbr.GetIdx())
                # Only need one valid candidate per –OH atom.
                break

    if not candidate_carbon_idxs:
        return False, "No valid aliphatic hydroxyl group found (or all -OH groups are in non-alcoholic contexts)"
    
    # Step 2: Build a graph of non-aromatic carbon atoms.
    # The graph is a dict mapping atom indices to a set of neighbor indices.
    carbon_graph = {}
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6 and not atom.GetIsAromatic():
            carbon_graph[atom.GetIdx()] = set()
            
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        if (a1.GetAtomicNum() == 6 and not a1.GetIsAromatic() and 
            a2.GetAtomicNum() == 6 and not a2.GetIsAromatic()):
            idx1, idx2 = a1.GetIdx(), a2.GetIdx()
            if idx1 in carbon_graph and idx2 in carbon_graph:
                carbon_graph[idx1].add(idx2)
                carbon_graph[idx2].add(idx1)
                
    if not carbon_graph:
        return False, "No non-aromatic carbon chain found."
    
    # BFS helper to compute distances.
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
        # Find farthest node and maximum distance.
        farthest_node = start
        max_distance = 0
        for node, dist in distances.items():
            if dist > max_distance:
                max_distance = dist
                farthest_node = node
        return farthest_node, distances
    
    # Compute the diameter (longest shortest–path) of a connected component.
    def component_diameter(component_nodes, graph):
        start = next(iter(component_nodes))
        node_a, _ = bfs(start, graph)
        node_b, distances = bfs(node_a, graph)
        return max(distances.values()) + 1  # chain length in atoms

    # Step 3: For each candidate carbon (attached to a valid –OH), get the connected component
    # and note its chain length (diameter).
    best_chain = 0
    processed_components = set()
    for candidate in candidate_carbon_idxs:
        if candidate not in carbon_graph:
            continue
        if candidate in processed_components:
            continue
        # Get the full connected component.
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
        if chain_length > best_chain:
            best_chain = chain_length

    # Step 4: Set required chain length threshold.
    # For small molecules (few heavy atoms) a chain as short as 3 carbons qualifies,
    # but for larger molecules we require at least 7 carbons.
    if mol.GetNumHeavyAtoms() > 20:
        min_chain = 7
    else:
        min_chain = 3

    if best_chain < min_chain:
        return False, f"Found a contiguous aliphatic carbon chain of {best_chain} atoms; need at least {min_chain}."

    # Step 5: Count the number of candidate hydroxyl groups detected in a similar manner.
    # (Too many may indicate a more functionalized, complex molecule.)
    candidate_OH_count = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 8:
            continue
        h_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 1]
        if len(h_neighbors) != 1:
            continue
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6 and not nbr.GetIsAromatic():
                is_carbonyl = False
                for cnbr in nbr.GetNeighbors():
                    if cnbr.GetAtomicNum() == 8 and cnbr.GetIdx() != atom.GetIdx():
                        bond = mol.GetBondBetweenAtoms(nbr.GetIdx(), cnbr.GetIdx())
                        if bond is not None and bond.GetBondType() == Chem.BondType.DOUBLE:
                            is_carbonyl = True
                            break
                if not is_carbonyl:
                    candidate_OH_count += 1
                    break
    if candidate_OH_count > 2:
        return False, f"Too many hydroxyl groups ({candidate_OH_count}) for a simple fatty alcohol."

    return True, (f"Found a contiguous aliphatic carbon chain of {best_chain} atoms "
                  f"attached to a valid hydroxyl group, with {candidate_OH_count} candidate -OH group(s).")

# For optional testing – remove this section in production.
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
        "C(CCCCCC)O",  # heptan-1-ol
    ]
    for s in test_smiles:
        result, reason = is_fatty_alcohol(s)
        print(f"SMILES: {s}\nResult: {result}\nReason: {reason}\n{'-'*60}")