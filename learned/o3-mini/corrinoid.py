"""
Classifies: CHEBI:33913 corrinoid
"""
"""
Classifies: Corrinoid
Definition: A derivative of the corrin nucleus, which contains four reduced or partly reduced
pyrrole rings joined in a macrocycle by three =C- groups and one direct carbon-carbon bond linking alpha positions.
This improved heuristic method first collects candidate five-membered rings having exactly one nitrogen.
Then, it builds a graph connecting rings that share an atom.
Finally, it attempts to “grow” a macrocycle from the connected candidate rings and verifies that the union of
atoms is itself part of a large ring (macrocycle) that is roughly the right size.
If any step fails, the molecule is rejected.
"""
from rdkit import Chem

def is_corrinoid(smiles: str):
    """
    Determines if a molecule is a corrinoid derivative based on its SMILES string.
    This function uses a heuristic approach: it looks for at least 4 five-membered rings (pyrrole-like) each having exactly one nitrogen, 
    then checks that they are all interconnected and form (or are part of) one macrocycle.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule is likely a corrinoid, False otherwise.
        str: Reason for classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()  # Each ring is a tuple of atom indices
    
    candidate_rings = []
    # Heuristic: Look for rings of size 5 (pyrrole-like) with exactly one nitrogen atom.
    for ring in atom_rings:
        if len(ring) == 5:
            n_count = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 7)
            if n_count == 1:
                candidate_rings.append(set(ring))
    
    if len(candidate_rings) < 4:
        return False, f"Found only {len(candidate_rings)} pyrrole-like five-membered rings (need at least 4)."
    
    # Build a graph connecting candidate rings that share at least one atom.
    # We require that the shared atom most likely is a bridging carbon or nitrogen (atomic number 6 or 7).
    n = len(candidate_rings)
    adj = {i: set() for i in range(n)}
    for i in range(n):
        for j in range(i+1, n):
            shared = candidate_rings[i].intersection(candidate_rings[j])
            if shared:
                # Check that at least one shared atom is C or N.
                if any(mol.GetAtomWithIdx(idx).GetAtomicNum() in (6,7) for idx in shared):
                    adj[i].add(j)
                    adj[j].add(i)
    
    # Find connected components in the graph of candidate rings.
    visited = set()
    comps = []
    for i in range(n):
        if i not in visited:
            stack = [i]
            comp = set()
            while stack:
                cur = stack.pop()
                if cur in comp:
                    continue
                comp.add(cur)
                for neigh in adj[cur]:
                    if neigh not in comp:
                        stack.append(neigh)
            visited |= comp
            comps.append(comp)
    
    # Now require at least one connected component that contains at least 4 candidate rings.
    good_comp = None
    for comp in comps:
        if len(comp) >= 4:
            good_comp = comp
            break
    if good_comp is None:
        return False, "The candidate pyrrole-like rings are not sufficiently interconnected into one macrocycle."
    
    # "Grow" the macrocycle: take the union of all atoms from these candidate rings.
    union_atoms = set()
    for i in good_comp:
        union_atoms |= candidate_rings[i]
    
    # To verify that the union forms part of one macrocycle, we use RDKit to see if these atoms
    # appear together in any ring in the molecule.
    union_atoms_list = list(union_atoms)
    submol = Chem.PathToSubmol(mol, union_atoms_list)
    sub_ring_info = submol.GetRingInfo().AtomRings()
    if not sub_ring_info:
        return False, "The union of candidate rings does not appear in any ring structure."
    # We expect one of the rings in the submolecule to be nearly as large as the union.
    max_ring_size = max(len(r) for r in sub_ring_info)
    # Allow a tolerance of 2 atoms (bridging atoms might not all be in candidate rings).
    if max_ring_size < len(union_atoms) - 2:
        return False, ("The candidate rings are not part of a single contiguous macrocycle; " +
                       f"union has {len(union_atoms)} atoms but largest ring in substructure has {max_ring_size}.")
    
    # In addition, to rule out extremely small systems, we can require that
    # the union macrocycle has a minimum size (e.g., at least 12 atoms).
    if len(union_atoms) < 12:
        return False, f"Macrocycle too small ({len(union_atoms)} atoms) to be a corrin nucleus."
    
    # If all conditions are met, we classify the molecule as a corrinoid.
    return True, "Molecule exhibits a corrin-nucleus-like pattern with 4 interconnected pyrrole-like rings forming a macrocycle."

# If the heuristic is too strained in any case, one might return (None, None) instead.
# (For example, an even more lenient approach might call return None, None if unexpected patterns are found.)