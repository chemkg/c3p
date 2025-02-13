"""
Classifies: CHEBI:47909 3-oxo-Delta(4) steroid
"""
"""
Classifies: 3-oxo-Delta(4) steroid – a 3-oxo steroid with an alpha,beta-unsaturated ketone (enone)
located on a fused steroid nucleus. The steroid nucleus is expected to be a contiguous fused ring system 
composed mostly of carbons and, in classical steroids, comprises three six‐membered rings and one five‐membered ring.
This implementation accepts candidate clusters that have at least three rings with at least one five‐membered ring 
and at least two six‐membered rings, and where at least 80% of the atoms are carbon.
"""

from rdkit import Chem

def is_3_oxo_Delta_4__steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-Delta(4) steroid based on its SMILES.
    The molecule must have a fused steroid nucleus (a fused ring system whose rings are mostly carbons and
    resemble that of a steroid – ideally containing at least one 5-membered ring and at least two 6-membered rings)
    as well as an enone motif (alpha,beta-unsaturated ketone) that is located on that nucleus.
    
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
    
    # Get ring information.
    ring_info = mol.GetRingInfo()
    all_rings = ring_info.AtomRings()
    if not all_rings:
        return False, "No rings found in the molecule"
    
    # Convert each ring (tuple of atom indices) to a set.
    rings = [set(r) for r in all_rings]
    n = len(rings)
    
    # Use a union-find algorithm to group rings that are fused.
    # Two rings are considered fused if they share at least 2 atoms.
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
    
    candidate_cluster = None
    candidate_reason = ""
    # Loop over each fused cluster.
    for group in groups.values():
        # We now relax the requirement: accept clusters with >= 3 rings.
        if len(group) >= 3:
            # Count ring sizes in this group.
            ring_sizes = [len(rings[i]) for i in group]
            count_5 = sum(1 for s in ring_sizes if s == 5)
            count_6 = sum(1 for s in ring_sizes if s == 6)
            # Require at least one five-membered and two six-membered rings.
            if count_5 >= 1 and count_6 >= 2:
                # Compute union of atom indices for this cluster.
                cluster_atoms = set()
                for ring_index in group:
                    cluster_atoms.update(rings[ring_index])
                # Check that the candidate nucleus consists mostly of carbon.
                n_carbons = sum(1 for idx in cluster_atoms if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
                if len(cluster_atoms) > 0 and (n_carbons / len(cluster_atoms)) >= 0.8:
                    candidate_cluster = cluster_atoms
                    break  # pick first acceptable cluster
                else:
                    candidate_reason = "Fused ring cluster found but does not consist mostly of carbons"
    
    if candidate_cluster is None:
        # If none satisfy the nucleus pattern.
        if candidate_reason:
            return False, candidate_reason
        return False, "The molecule does not have a fused ring cluster with a steroid nucleus pattern (>=3 fused rings with at least 1 five-membered and 2 six-membered rings, mostly carbon)"
    
    # Define a SMARTS pattern for the enone motif (alpha,beta-unsaturated ketone).
    # The pattern here requires a ring carbon (in a ring) bearing a carbonyl group,
    # connected to a carbon-carbon double bond (with both atoms in rings).
    enone_pattern = Chem.MolFromSmarts("[#6;R](=O)[#6;R]=[#6;R]")
    if enone_pattern is None:
        return False, "Error in generating SMARTS pattern for enone motif"
    
    enone_matches = mol.GetSubstructMatches(enone_pattern)
    if not enone_matches:
        return False, "No alpha,beta-unsaturated ketone (enone) motif found in the molecule"
    
    # Check that at least one enone motif is within the candidate steroid nucleus.
    for match in enone_matches:
        # We check that the carbonyl atom (first atom in the match) is in the candidate nucleus.
        if match[0] in candidate_cluster:
            return True, "Contains a fused steroid nucleus (mostly carbons, with >=1 five-membered and >=2 six-membered rings) and an enone motif on it (3-oxo-Delta(4) steroid)"
    
    return False, "An enone motif was found but not located within the fused steroid nucleus"