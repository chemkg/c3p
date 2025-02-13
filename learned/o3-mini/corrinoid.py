"""
Classifies: CHEBI:33913 corrinoid
"""
"""
Classifies: Corrinoid
Definition: A derivative of the corrin nucleus, which contains four reduced or partly reduced 
pyrrole rings joined in a macrocycle by three =C- groups and one direct carbon-carbon bond linking alpha positions.

Heuristic approach in this code:
 1. Parse the molecule from its SMILES.
 2. Identify candidate five-membered rings having exactly one nitrogen (which mimics a pyrrole ring).
 3. Take the union of all atoms in these candidate rings.
 4. Build a connectivity graph among these atoms (using bonds that connect two atoms in the union) and ensure it is connected.
 5. Check that the molecule contains a cobalt atom (atomic number 27) which is typical in corrinoids.
 6. Require that there are at least 4 candidate rings and that the union of their atoms is sufficiently large.
If any step fails, the function returns False along with a reason.
"""

from rdkit import Chem

def is_corrinoid(smiles: str):
    """
    Determines if a molecule is a corrinoid derivative based on its SMILES string using a heuristic approach.
    
    The technique is as follows:
      - It gathers candidate five-membered rings (as provided by RDKit ring info) that contain exactly one nitrogen atom,
        which is typical of pyrrole rings.
      - It then takes the union of the atoms in these candidate rings and constructs a connectivity graph (by considering bonds
        that connect two atoms in that union). If the union is not connected, we reject the molecule.
      - Additionally, because most corrinoids contain a cobalt (Co) ion (atomic number 27), we require its presence.
      - Some minimal size (number of atoms in the union) is also enforced.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is likely a corrinoid, False otherwise.
        str: A reason for the resulting classification.
    """
    # Parse the SMILES string to a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # Retrieve ring information from the molecule.
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()  # Each ring is a tuple of atom indices.
    
    # Gather candidate rings that are five-membered and have exactly one nitrogen atom.
    candidate_rings = []
    for ring in atom_rings:
        if len(ring) == 5:
            n_count = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 7)
            # In a five-membered ring, exactly 1 nitrogen is typical of a pyrrole ring.
            if n_count == 1:
                candidate_rings.append(set(ring))
    
    if len(candidate_rings) < 4:
        return False, f"Found only {len(candidate_rings)} candidate five-membered pyrrole-like rings (need at least 4)."
    
    # Build the union of all atoms in candidate rings.
    union_atoms = set().union(*candidate_rings)
    
    # Build a connectivity graph among the union atoms (nodes are atom indices; an edge exists if both atoms share a bond).
    graph = {idx: set() for idx in union_atoms}
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtomIdx()
        a2 = bond.GetEndAtomIdx()
        if a1 in union_atoms and a2 in union_atoms:
            graph[a1].add(a2)
            graph[a2].add(a1)
    
    # Check connectivity: ensure that the subgraph induced by union_atoms is connected.
    visited = set()
    def dfs(node):
        stack = [node]
        while stack:
            current = stack.pop()
            if current not in visited:
                visited.add(current)
                stack.extend(graph[current] - visited)
    # Start from one arbitrary atom in the union.
    start = next(iter(union_atoms))
    dfs(start)
    if visited != union_atoms:
        return False, "Candidate rings are not connected into one contiguous macrocycle."
    
    # Check for the presence of a cobalt atom (atomic number 27) which is a hallmark in most corrinoids.
    has_cobalt = any(atom.GetAtomicNum() == 27 for atom in mol.GetAtoms())
    if not has_cobalt:
        return False, "No cobalt atom found in the molecule; it may not be a corrinoid."
    
    # Optional: require that the union of candidate ring atoms is of a minimum size.
    if len(union_atoms) < 12:
        return False, f"Union of candidate rings is too small ({len(union_atoms)} atoms) to represent a corrin macrocycle."
    
    return True, "Molecule exhibits a corrin-nucleus-like pattern with four interconnected pyrrole-like rings and a cobalt center."

# Note: The heuristic method here is simple and may need further refinement for real world applications.