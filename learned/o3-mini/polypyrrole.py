"""
Classifies: CHEBI:38077 polypyrrole
"""
#!/usr/bin/env python
"""
Classifies: A compound composed of two or more pyrrole units.
A pyrrole unit here is defined as a 5‐membered aromatic ring 
with exactly one nitrogen and four carbon atoms.
We first extract rings from the molecule (via GetRingInfo), confirm which rings 
meet the pyrrole criteria, and then group nearby rings into connected clusters.
Two rings are considered “connected” if the ring plus its immediate neighbors 
(over its bonds) touches the other ring. In addition, if two rings are directly fused 
(i.e. share at least one atom) we lower the fraction threshold.
Finally, for the largest connected cluster found we compute 
the fraction f = (# unique heavy atoms in the pyrrole units)/(# heavy atoms in molecule).
If the cluster has:
  • 2 pyrrole units: if they are fused, we require f >= 0.20, otherwise we require f >= 0.35.
  • ≥3 pyrrole units: we require f >= 0.40.
If at least one cluster meets the above criteria the molecule qualifies as a polypyrrole.
Note: This rule‐based method remains heuristic and may miss some examples or yield false positives.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_polypyrrole(smiles: str):
    """
    Determines if a molecule is a polypyrrole – i.e. it contains at least two pyrrole units 
    (as defined: 5-membered aromatic rings with one nitrogen and four carbons) arranged in a connected network.
    
    The algorithm proceeds as follows:
      1. Parse and sanitize the molecule.
      2. Count heavy atoms (non‐H) for later fraction computation.
      3. Identify all rings (from GetRingInfo) and keep only those 5-membered, aromatic rings with exactly one N.
      4. For each candidate pyrrole ring, expand its atom set to include all directly-bonded neighbors.
      5. Cluster pyrrole rings as “connected” if one ring’s expanded set overlaps the other ring’s atom set.
      6. Mark clusters as “fused” if any two rings in the cluster share at least one atom.
      7. For each cluster, compute the union of the actual ring atoms and its heavy atom fraction.
         • For a two‐ring cluster: require f >= 0.20 if fused, else f >= 0.35.
         • For clusters of 3 or more rings, require f >= 0.40.
      8. Return True if at least one cluster meets the criteria.
      
    Args:
       smiles (str): SMILES string of the molecule.
       
    Returns:
       bool: True if the molecule qualifies as a polypyrrole, False otherwise.
       str: Detailed reason for the decision.
    """
    # Parse SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # Sanitize the molecule to assign aromaticity etc.
    try:
        Chem.SanitizeMol(mol)
    except Exception as e:
        return False, f"Sanitization failed: {e}"
    
    # Count heavy atoms (non-hydrogens) for fraction computation.
    heavy_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1]
    if not heavy_atoms:
        return False, "No heavy atoms in molecule."
    total_heavy = len(heavy_atoms)
    
    # Extract rings from the molecule.
    ring_info = mol.GetRingInfo()
    candidate_rings = []
    for ring in ring_info.AtomRings():
        # Only 5-membered rings are candidate for a pyrrole unit.
        if len(ring) != 5:
            continue
        # Check that every ring atom is aromatic.
        if not all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            continue
        # Count nitrogen atoms.
        n_count = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 7)
        if n_count == 1:
            candidate_rings.append(set(ring))
    
    if len(candidate_rings) < 2:
        return False, f"Found {len(candidate_rings)} pyrrole ring(s); at least 2 are required."
    
    # For connectivity we expand each ring's atom set by including neighbors.
    expanded_rings = []
    for ring in candidate_rings:
        expanded = set(ring)
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            for nbr in atom.GetNeighbors():
                expanded.add(nbr.GetIdx())
        expanded_rings.append(expanded)
    
    # Cluster the rings via union-find; two rings are “connected” if the expanded set of one overlaps
    # the actual atom set of the other.
    n_rings = len(candidate_rings)
    parent = list(range(n_rings))
    
    def find(i):
        # Find with path compression.
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i
    
    def union(i, j):
        ri, rj = find(i), find(j)
        if ri != rj:
            parent[rj] = ri
    
    for i in range(n_rings):
        for j in range(i+1, n_rings):
            if expanded_rings[i] & candidate_rings[j] or expanded_rings[j] & candidate_rings[i]:
                union(i, j)
                
    # Group rings by their cluster roots.
    clusters = {}
    for i in range(n_rings):
        root = find(i)
        clusters.setdefault(root, []).append(i)
    
    # Evaluate each cluster.
    # For each cluster, compute:
    #   a) The number of pyrrole rings (n_ring)
    #   b) The union of the ring atoms (actual atoms in the candidate rings)
    #   c) Whether the cluster is "fused" (any two rings share at least one atom directly)
    cluster_found = False
    for cluster_indices in clusters.values():
        n_ring = len(cluster_indices)
        if n_ring < 2:
            continue  # we need at least 2 pyrrole units in a cluster
        union_atoms = set()
        rings_in_cluster = []
        for idx in cluster_indices:
            ring_atoms = candidate_rings[idx]
            union_atoms |= ring_atoms
            rings_in_cluster.append(ring_atoms)
        fraction = len(union_atoms) / total_heavy
        
        # Determine if any pair is directly fused (i.e. share at least one atom)
        fused = False
        for i in range(len(rings_in_cluster)):
            for j in range(i+1, len(rings_in_cluster)):
                if rings_in_cluster[i] & rings_in_cluster[j]:
                    fused = True
                    break
            if fused:
                break
                
        # Set threshold based on cluster type.
        if n_ring == 2:
            threshold = 0.20 if fused else 0.35
        else:
            threshold = 0.40

        # If cluster meets threshold, we can classify as polypyrrole.
        if fraction >= threshold:
            reason = (f"Connected cluster has {n_ring} pyrrole unit(s) "
                      f"(fused: {fused}) covering {fraction:.2f} of heavy atoms "
                      f"(threshold {threshold}). Qualifies as a polypyrrole.")
            return True, reason
            
    # If no cluster qualifies, return details based on the best (largest or highest fraction) cluster.
    best_cluster = None
    best_n = 0
    best_frac = 0.0
    for cluster_indices in clusters.values():
        n_ring = len(cluster_indices)
        union_atoms = set()
        for idx in cluster_indices:
            union_atoms |= candidate_rings[idx]
        frac = len(union_atoms) / total_heavy
        if n_ring > best_n or (n_ring == best_n and frac > best_frac):
            best_n = n_ring
            best_frac = frac
            best_cluster = cluster_indices
    return False, (f"No connected pyrrole cluster met the threshold. "
                   f"Best cluster had {best_n} pyrrole unit(s) covering {best_frac:.2f} "
                   f"of heavy atoms.")

# Example usage for testing:
if __name__ == "__main__":
    test_cases = [
        # True positive examples:
        ("CN(C)CCNC(=O)CCC1=[N+]2C(C=C1)=Cc1ccc(-c3ccc[nH]3)n1[B-]2(F)F", "DND-99 dye"),
        ("O(C=1C(N=C(C1)C=2NC=CC2)=CC=3NC(C=4NC=CC4)=CC3OC)C", "Purple pigment"),
        ("c1cc2cc3ccc([nH]3)c3ccc(cc4ccc(cc5ccc(cc1n2)[nH]5)n4)[nH]3", "Sapphyrin"),
        ("[H]C(=C([H])C1=[N+]2C(C=C1)=Cc1ccc(-c3ccc[nH]3)n1[B-]2(F)F)c1ccc(OCC(=O)NCCCCCC(=O)ON2C(=O)CCC2=O)", "BODIPY 650/665-X"),
        ("O(C1=C2NC(=C1)C=3NC(CCCCCCCCCC4=NC(=C2)C=C4)=CC3)C", "Nonylprodigiosin"),
        
        # False positive example:
        ("ClC1=CC2=C(N([C@@H]3O[C@@H]([C@@H](O)[C@@H]([C@H]3O)O)C)C4=C2C=CC5=C4NC=6C=CC=CC56)C=C1", "Tjipanazole C1"),
        
        # False negative examples (should be rejected because fewer than 2 connected pyrrole rings detected):
        ("C=1(C2=C(NC1CC=3NC(C(=O)[H])=C(C)C3CC)/C(=C/4\\N=C(C[C@@]5(C(C)=C(C=C)C(=O)N5)[H])[C@@H](C)[C@@H]4CCC(=O)O)/[C@@H](C(=O)OC)C2=O", 
         "(1R)-primary fluorescent chlorophyll catabolite"),
    ]
    
    for smi, name in test_cases:
        result, reason = is_polypyrrole(smi)
        print(f"Name: {name}\nResult: {result}\nReason: {reason}\n{'-'*60}")