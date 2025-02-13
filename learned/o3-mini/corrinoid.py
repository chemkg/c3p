"""
Classifies: CHEBI:33913 corrinoid
"""
"""
Classifies: Corrinoid
Definition: A derivative of the corrin nucleus, which contains four reduced or partly reduced 
pyrrole rings joined in a macrocycle by three =C- groups and one direct carbon-carbon bond linking 
alpha positions.

Heuristic approach in this code:
 1. Parse the molecule from its SMILES.
 2. Identify candidate rings as five-membered rings that have exactly one nitrogen atom (as expected for a pyrrole ring).
 3. Instead of taking the union of atoms from each candidate ring and checking connectivity, we build 
    a connectivity graph where each candidate ring is a node and an edge exists between two rings if any 
    atom in one ring is directly bonded (in the full molecule) to any atom of the other ring.
 4. We then check that the candidate ring graph is connected (all nodes can be reached from one another).
 5. Also require that there are at least 4 candidate rings.
 6. Finally, ensure that the molecule contains a cobalt atom (Co, atomic number 27) which is a hallmark in corrinoids.
"""

from rdkit import Chem

def is_corrinoid(smiles: str):
    """
    Determines if a molecule is a corrinoid derivative based on its SMILES string using a heuristic approach.
    
    The technique is as follows:
      - Identify candidate five-membered rings (from RDKit ring info) that contain exactly one nitrogen atom.
      - Build a graph whose nodes represent these candidate rings.
        Two candidate rings are connected by an edge if at least one atom in one ring is directly bonded
        to any atom in the other ring.
      - Check that there are at least four candidate rings and that these rings form one connected network.
      - Ensure that a cobalt atom (atomic number 27) is present.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is likely a corrinoid, False otherwise.
        str: A reason for the resulting classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # Retrieve ring information
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()  # Each ring is a tuple of atom indices.
    
    # Gather candidate rings: five-membered rings with exactly one nitrogen atom
    candidate_rings = []
    for ring in atom_rings:
        if len(ring) == 5:
            n_count = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 7)
            # For a pyrrole, exactly one nitrogen is expected in the five-membered ring.
            if n_count == 1:
                candidate_rings.append(set(ring))
    
    if len(candidate_rings) < 4:
        return False, f"Found only {len(candidate_rings)} candidate five-membered pyrrole-like rings (need at least 4)."
    
    # Build a graph where nodes correspond to candidate rings (indexed by their position in candidate_rings list)
    # Two candidate rings are connected if any atom in one ring is bonded to any atom in the other (using the full molecule bonds).
    ring_graph = {i: set() for i in range(len(candidate_rings))}
    # Iterate over all pairs of candidate rings
    for i in range(len(candidate_rings)):
        for j in range(i+1, len(candidate_rings)):
            # Check all pairs of atoms between ring i and ring j
            connected = False
            for a in candidate_rings[i]:
                for b in candidate_rings[j]:
                    if mol.GetBondBetweenAtoms(a, b) is not None:
                        connected = True
                        break
                if connected:
                    break
            if connected:
                ring_graph[i].add(j)
                ring_graph[j].add(i)
    
    # Check connectivity in the ring graph (using DFS)
    visited = set()
    def dfs(node):
        stack = [node]
        while stack:
            current = stack.pop()
            if current not in visited:
                visited.add(current)
                stack.extend(ring_graph[current] - visited)
    dfs(0)
    if len(visited) != len(candidate_rings):
        return False, "Candidate rings are not connected into one contiguous macrocycle."
    
    # Check for the presence of a cobalt (Co) atom (atomic number 27). Corrinoids usually contain a cobalt center.
    has_cobalt = any(atom.GetAtomicNum() == 27 for atom in mol.GetAtoms())
    if not has_cobalt:
        return False, "No cobalt atom found in the molecule; it may not be a corrinoid."
    
    return True, "Molecule exhibits a corrinoid signature with at least 4 interconnected pyrrole-like rings and a cobalt center."