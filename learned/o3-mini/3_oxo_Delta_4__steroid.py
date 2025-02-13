"""
Classifies: CHEBI:47909 3-oxo-Delta(4) steroid
"""
"""
Classifies: 3-oxo-Delta(4) steroid – a 3-oxo steroid with an alpha,beta-unsaturated ketone (enone)
located on a fused steroid nucleus. The steroid nucleus must be built of a fused system with at least 4 rings,
with a characteristic ring-size distribution (at least one five-membered ring and at least three six-membered rings).
"""
from rdkit import Chem

def is_3_oxo_Delta_4__steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-Delta(4) steroid based on its SMILES.
    The molecule must have a steroid nucleus (fused rings with at least one 5-membered and three 6-membered rings)
    as well as an enone motif (alpha,beta-unsaturated ketone) located on that nucleus.
    
    Args:
        smiles (str): SMILES string for the molecule.
        
    Returns:
        bool: True if the criteria are met, False otherwise.
        str: A reason for the classification result.
    """
    # Parse and sanitize the molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    try:
        Chem.SanitizeMol(mol)
    except Exception as e:
        return False, f"Sanitization error: {str(e)}"
        
    # Get ring information from the molecule.
    ring_info = mol.GetRingInfo()
    all_rings = ring_info.AtomRings()
    if not all_rings:
        return False, "No rings found in the molecule"
    
    # Convert each ring (tuple of atom indices) to a set for easier intersection handling.
    rings = [set(r) for r in all_rings]
    
    # Use a union-find algorithm to group rings that are fused.
    # Two rings are fused if they share at least 2 atoms.
    n = len(rings)
    parent = list(range(n))
    
    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i
    
    def union(i, j):
        pi = find(i)
        pj = find(j)
        if pi != pj:
            parent[pj] = pi
    
    for i in range(n):
        for j in range(i + 1, n):
            if len(rings[i].intersection(rings[j])) >= 2:
                union(i, j)
    
    # Group rings by their connected component.
    groups = {}
    for i in range(n):
        root = find(i)
        groups.setdefault(root, []).append(i)
    
    # Look for a fused cluster (group) that is likely to represent a steroid nucleus.
    # For a typical steroid nucleus, we expect at least four rings, with three six-membered rings and one five-membered ring.
    candidate_cluster = None
    for group in groups.values():
        if len(group) >= 4:
            # For the rings in this group, count the ring sizes.
            ring_sizes = [len(rings[i]) for i in group]
            count_5 = sum(1 for s in ring_sizes if s == 5)
            count_6 = sum(1 for s in ring_sizes if s == 6)
            if count_5 >= 1 and count_6 >= 3:
                candidate_cluster = group
                break
    if candidate_cluster is None:
        return False, "The molecule does not have a fused ring cluster with the steroid nucleus pattern (at least one five-membered and three six-membered rings)"
    
    # Compute the union of atom indices from the candidate steroid cluster.
    steroid_atoms = set()
    for ring_index in candidate_cluster:
        steroid_atoms.update(rings[ring_index])
    
    # Define a SMARTS pattern for the enone motif (alpha,beta-unsaturated ketone).
    # The pattern looks for a ring carbon (in a ring) with a carbonyl (C(=O)) 
    # attached to another ring carbon that is double-bonded to a third ring carbon.
    enone_pattern = Chem.MolFromSmarts("[#6;R](=O)[#6;R]=[#6;R]")
    if enone_pattern is None:
        return False, "Error in generating SMARTS pattern for enone motif"
    
    enone_matches = mol.GetSubstructMatches(enone_pattern)
    if not enone_matches:
        return False, "No alpha,beta-unsaturated ketone (enone) motif found in the molecule"
    
    # Require that at least one enone is located within the candidate steroid nucleus.
    # We use the location of the carbonyl atom (the first atom in the match).
    for match in enone_matches:
        carbonyl_atom_idx = match[0]
        if carbonyl_atom_idx in steroid_atoms:
            return True, "Contains a fused steroid nucleus (with appropriate ring sizes) and an enone motif in it (3-oxo-Delta(4) steroid)"
    
    return False, "An enone motif was found but not located within the fused steroid nucleus"

# Example (uncomment below to run tests):
# test_smiles = [
#     # True examples:
#     "[H][C@@]12C[C@H](C)[C@](O)(C(=O)CO)[C@@]1(C)C[C@H](O)[C@@]1(F)[C@@]2([H])CCC2=CC(=O)C=C[C@]12C",  # betamethasone
#     "C1C(C=C2[C@](C1)([C@@]3([C@@](CC2)([C@@]4([H])[C@@](CC3)(C)[C@H]([C@@H](C4)O)O)[H])[H])C)=O",  # 16alpha-hydroxytestosterone
#     # False positive example (non-steroid with enone motif)
#     "C[C@H]1O[C@H](C[C@@H](O)[C@@H]1O)c1ccc2C(=O)C3=C([C@H](O)C[C@]4(O)C[C@@](C)(O)CC(=O)[C@]34O)C(=O)c2c1O"
# ]
# for s in test_smiles:
#     flag, reason = is_3_oxo_Delta_4__steroid(s)
#     print(flag, reason)