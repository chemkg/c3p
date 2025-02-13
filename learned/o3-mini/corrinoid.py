"""
Classifies: CHEBI:33913 corrinoid
"""
"""
Classifies: Corrinoid
Definition: A derivative of the corrin nucleus, which contains four reduced or partly reduced
pyrrole rings joined in a macrocycle by three =C- groups and one direct carbon-carbon bond linking alpha positions.
This heuristic function seeks to identify four five-membered rings (pyrrole-like rings with one nitrogen)
that are connected.
"""
from rdkit import Chem

def is_corrinoid(smiles: str):
    """
    Determines if a molecule is a corrinoid derivative based on its SMILES string.
    The approach heuristically looks for four five-membered rings that each contain exactly one nitrogen
    (as expected for pyrrole rings in the corrin nucleus) and are connected as part of a macrocycle.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule is likely a corrinoid, False otherwise.
        str: Reason for classification.
    """
    # Parse SMILES string (if invalid, return immediately)
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Get ring information from the molecule
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()  # This returns tuples of atom indices for each ring

    candidate_rings = []
    # Heuristically select rings that are likely pyrrole-like:
    # Expecting rings of size 5 with exactly one nitrogen atom.
    for ring in atom_rings:
        if len(ring) == 5:
            n_count = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 7)
            if n_count == 1:
                candidate_rings.append(set(ring))
    
    if len(candidate_rings) < 4:
        return False, f"Found only {len(candidate_rings)} pyrrole-like five-membered rings (need at least 4)."
    
    # Next, check interconnectivity of these candidate rings.
    # We build a simple graph where each candidate ring is a node, and an edge exists if two rings share at least one atom.
    n = len(candidate_rings)
    adj = {i: set() for i in range(n)}
    for i in range(n):
        for j in range(i+1, n):
            if candidate_rings[i].intersection(candidate_rings[j]):
                # There is a shared atom (which ideally is a carbon bridging the rings)
                # We check that the shared atom is carbon (atomic number 6) if possible.
                shared = candidate_rings[i].intersection(candidate_rings[j])
                has_carbon = any(mol.GetAtomWithIdx(idx).GetAtomicNum() == 6 for idx in shared)
                if has_carbon:
                    adj[i].add(j)
                    adj[j].add(i)
    
    # To be part of one macrocycle, the four rings should be connected.
    # We try to find a connected component of candidate rings of size at least 4.
    visited = set()
    def dfs(node):
        stack = [node]
        comp = set()
        while stack:
            cur = stack.pop()
            if cur in comp:
                continue
            comp.add(cur)
            for neigh in adj[cur]:
                if neigh not in comp:
                    stack.append(neigh)
        return comp

    # Check all components.
    connected = False
    for i in range(n):
        comp = dfs(i)
        if len(comp) >= 4:
            connected = True
            break

    if not connected:
        return False, "The candidate pyrrole-like rings are not sufficiently interconnected to form a corrin nucleus."
    
    # The molecule has at least 4 connected pyrrole-like rings.
    return True, "Molecule exhibits a corrin-nucleus-like pattern with 4 interconnected pyrrole-like rings."