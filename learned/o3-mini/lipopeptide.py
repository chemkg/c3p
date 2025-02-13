"""
Classifies: CHEBI:46895 lipopeptide
"""
#!/usr/bin/env python3
"""
Classifies: Lipopeptide
Definition: A compound consisting of a peptide (a molecule containing at least one amide bond)
            with an attached lipid. The lipid portion is heuristically defined as having at least 
            one contiguous acyclic chain of non-aromatic carbon atoms of a minimum length (here 6 atoms).
            
Below is our improved algorithm:
  (1) Verify that the molecule contains at least one peptide (amide) bond (using the SMARTS "C(=O)N").
  (2) Look for “lipid tails” by first extracting the set of atoms that are:
         - Carbon atoms (atomic number 6)
         - Not aromatic and not in a ring (i.e. not part of cyclic or aromatic moieties)
      Then we define the subgraph induced by these atoms and compute the longest (graph) path within 
      each connected component. If any such longest path is at least 6 atoms long, we flag the molecule 
      as having a lipid tail.
      
This revised method is tuned so that some false positives (where a peptide simply has a chain of carbons
internally) are avoided while not missing known lipopeptides such as colistimethate B.
"""

from rdkit import Chem
from collections import deque

def _bfs_longest_path(graph, start):
    """Perform a BFS from start node in an undirected graph (dict) to find the farthest node and its distance (in number of edges)."""
    visited = {start: 0}
    queue = deque([start])
    farthest_node = start
    max_dist = 0
    while queue:
        current = queue.popleft()
        dist = visited[current]
        if dist > max_dist:
            max_dist = dist
            farthest_node = current
        for neighbor in graph.get(current, []):
            if neighbor not in visited:
                visited[neighbor] = dist + 1
                queue.append(neighbor)
    return farthest_node, max_dist

def _diameter_of_component(component, graph):
    """
    Given a list of node indices (component) and the subgraph (dict mapping node:index to neighbor list),
    compute the diameter (i.e. the length in number of nodes of the longest path).
    (We perform two BFS passes: first from any node to get a far vertex, then from that vertex.)
    """
    if not component:
        return 0
    start = component[0]
    far_node, _ = _bfs_longest_path(graph, start)
    _, longest_edge_path = _bfs_longest_path(graph, far_node)
    # The number of nodes along the path is (longest_edge_path + 1)
    return longest_edge_path + 1

def is_lipopeptide(smiles: str):
    """
    Determines if a molecule is a lipopeptide based on its SMILES string.
    A lipopeptide is defined as a compound that contains a peptide portion
    (indicated by at least one amide (C(=O)N) bond) and an attached lipid tail.
    
    The lipid tail is approximated by the presence of at least one contiguous linear chain 
    of non-aromatic, acyclic carbon atoms of a minimum length (here set to 6 atoms).
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a lipopeptide, False otherwise.
        str: Reason for the classification.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # (1) Check for peptide (amide) bond: use the SMARTS pattern "C(=O)N"
    amide_smarts = "C(=O)N"
    amide_pattern = Chem.MolFromSmarts(amide_smarts)
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "No amide (peptide) bonds found in the molecule."
    
    # (2) Look for a lipid chain.
    # Heuristic: We search for a contiguous chain (subgraph) of carbon atoms that are
    #         non-aromatic and not in any ring.
    #         We only follow C–C bonds.
    #         We require that the longest linear chain (diameter of the graph) is at least a threshold.
    lipid_chain_threshold = 6  # relaxed threshold from 8 to 6 atoms to reduce false negatives
    # Build a set of atoms eligible for a lipid chain: carbon atoms (atomic num == 6) that are
    # non-aromatic and not in a ring.
    eligible_atoms = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6 and (not atom.GetIsAromatic()) and (not atom.IsInRing()):
            eligible_atoms.append(atom.GetIdx())
    
    if not eligible_atoms:
        return False, "No eligible aliphatic (non-cyclic, non-aromatic) carbon atoms found for lipid tail."
    
    # Build an undirected graph among the eligible atoms, considering only bonds between two eligible atoms.
    graph = {}
    eligible_set = set(eligible_atoms)
    for atom_idx in eligible_atoms:
        atom = mol.GetAtomWithIdx(atom_idx)
        neighbors = []
        for neighbor in atom.GetNeighbors():
            if neighbor.GetIdx() in eligible_set:
                # Only follow bonds where both atoms are eligible (both carbons meeting our criteria)
                neighbors.append(neighbor.GetIdx())
        graph[atom_idx] = neighbors
    
    # Find connected components in the subgraph.
    visited = set()
    lipid_tail_found = False
    for atom_idx in eligible_atoms:
        if atom_idx in visited:
            continue
        # Do a DFS/BFS to gather the component
        component = []
        stack = [atom_idx]
        while stack:
            current = stack.pop()
            if current in visited:
                continue
            visited.add(current)
            component.append(current)
            for neighbor in graph.get(current, []):
                if neighbor not in visited:
                    stack.append(neighbor)
        # For this component, compute its diameter (longest chain in number of atoms)
        diameter = _diameter_of_component(component, graph)
        if diameter >= lipid_chain_threshold:
            lipid_tail_found = True
            break  # we only need one substantial lipid chain
    
    if not lipid_tail_found:
        return False, f"No contiguous chain found with at least {lipid_chain_threshold} carbon atoms (lipid tail)."
    
    return True, "Molecule contains peptide bonds and has a lipid tail (a contiguous chain of sufficient carbon atoms)."

# Example usage (for testing purposes)
if __name__ == '__main__':
    # Test one of the known lipopeptide structures (for example, surfactin A)
    surfactin_A = "[H][C@@]1(CCCCCCCC(C)C)CC(=O)N[C@@H](CCC(O)=O)C(=O)N[C@@H](CC(C)C)C(=O)N[C@H](CC(C)C)C(=O)N[C@@H](C(C)C)C(=O)N[C@@H](CC(O)=O)C(=O)N[C@H](CC(C)C)C(=O)N[C@@H](CC(C)C)C(=O)O1"
    result, reason = is_lipopeptide(surfactin_A)
    print("Surfactin A classification:", result, "-", reason)