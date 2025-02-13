"""
Classifies: CHEBI:51963 hopanoid
"""
"""
Classifies: A hopanoid, defined as a triterpenoid based on a hopane skeleton.
A hopanoid is (roughly) expected to have a saturated, fused polycyclic core 
with ≥5 fused rings (with at least one 5‐membered ring and three 6‐membered rings)
and a fused core that is around the typical triterpenoid carbon count (roughly 20–40 carbons)
with a very high sp3 fraction (i.e. nearly all carbons are saturated).
Note: Decorations (e.g. sugars, amino substituents) may lie outside the core.
If the classification is too uncertain, the function may return (None, None).
"""
from rdkit import Chem

def is_hopanoid(smiles: str):
    """
    Determines if a molecule belongs to the hopanoid class based on its SMILES string.
    
    Instead of relying on the overall carbon count (which may be skewed by non-core decorations),
    we first extract the largest fused ring cluster (considered the core) by checking rings that
    share ≥2 atoms. We then check that the cluster has at least 5 rings in total, contains at least
    one 5-membered ring and at least three 6-membered rings, and that the core (its union of ring atoms)
    has a carbon count in the expected range and an extremely high fraction of sp3 carbons.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a hopanoid, False otherwise.
        str: Explanation for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get ring information from RDKit (each ring is returned as a tuple of atom indices).
    ring_info = mol.GetRingInfo().AtomRings()
    n_rings = len(ring_info)
    if n_rings < 5:
        return False, f"Only {n_rings} rings detected; hopanoid core requires at least 5 fused rings."
    
    # Convert each ring into a set so we can test for fusions.
    ring_sets = [set(r) for r in ring_info]
    n = len(ring_info)
    # Build an adjacency list for rings sharing at least 2 atoms.
    adj = {i: set() for i in range(n)}
    for i in range(n):
        for j in range(i + 1, n):
            if len(ring_sets[i].intersection(ring_sets[j])) >= 2:
                adj[i].add(j)
                adj[j].add(i)
    
    # Find connected components (clusters of fused rings) using depth-first search.
    visited = set()
    clusters = []
    for i in range(n):
        if i in visited:
            continue
        stack = [i]
        comp = set()
        while stack:
            current = stack.pop()
            if current in comp:
                continue
            comp.add(current)
            for neighbor in adj[current]:
                if neighbor not in comp:
                    stack.append(neighbor)
        visited |= comp
        clusters.append(comp)
    
    # Use the largest fused ring cluster for further classification.
    largest_cluster = max(clusters, key=lambda comp: len(comp))
    if len(largest_cluster) < 5:
        return False, (f"Largest fused ring cluster contains only {len(largest_cluster)} rings; "
                       "at least 5 are required for a hopanoid core")
    
    # Check that the fused cluster contains the required ring sizes.
    # Count number of 5-membered and 6-membered rings in the cluster.
    cluster_ring_sizes = [len(ring_info[i]) for i in largest_cluster]
    num_5 = sum(1 for size in cluster_ring_sizes if size == 5)
    num_6 = sum(1 for size in cluster_ring_sizes if size == 6)
    if num_5 < 1:
        return False, "Fused ring cluster does not contain any 5-membered ring; at least one is required."
    if num_6 < 3:
        return False, f"Fused ring cluster contains only {num_6} six-membered rings; at least 3 are required."
    
    # Unite all atoms that are in any ring within the largest fused cluster. This is our 'core'.
    fused_atoms = set()
    for i in largest_cluster:
        fused_atoms |= ring_sets[i]
    
    # Count the carbon atoms within the core.
    core_carbons = [mol.GetAtomWithIdx(idx) for idx in fused_atoms if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6]
    core_carbon_count = len(core_carbons)
    if core_carbon_count < 20 or core_carbon_count > 40:
        return False, (f"Fused core contains {core_carbon_count} carbon atoms; "
                       "expected roughly 20-40 carbons for a hopane skeleton.")
    
    # Ensure that none of the atoms in the core are aromatic
    for idx in fused_atoms:
        if mol.GetAtomWithIdx(idx).GetIsAromatic():
            return False, "Fused ring core contains aromatic atoms; hopanoid cores are typically fully alicyclic (non‐aromatic)."
    
    # Calculate sp3 fraction among the core carbon atoms.
    n_sp3 = sum(1 for atom in core_carbons if atom.GetHybridization() == Chem.rdchem.HybridizationType.SP3)
    sp3_fraction_core = n_sp3 / core_carbon_count if core_carbon_count > 0 else 0
    if sp3_fraction_core < 0.85:
        return False, (f"Low sp3 fraction in core (fraction = {sp3_fraction_core:.2f}); "
                       "hopanoid cores are expected to be highly saturated.")
    
    # If reached here, the fused core meets many of the expected criteria.
    reason = (f"Contains a fused polycyclic core with {len(largest_cluster)} rings (ring sizes: {cluster_ring_sizes}), "
              f"core carbon count {core_carbon_count} and high core sp3 fraction ({sp3_fraction_core:.2f}), "
              "consistent with a hopanoid (hopane-type) skeleton.")
    return True, reason

# Example usage (for testing):
if __name__ == "__main__":
    # Known hopanoid example: (32R,33R,34R)-bacteriohopanetetrol
    test_smiles = "O[C@@H]([C@H](O)CO)[C@H](O)CC[C@H]([C@@H]1[C@@H]2[C@]([C@@H]3[C@@]([C@]4([C@@H]([C@@]5([C@H](C(CCC5)(C)C)CC4)C)CC3)C)(C)CC2)(C)CC1)C"
    status, explanation = is_hopanoid(test_smiles)
    print("Is hopanoid?", status)
    print("Reason:", explanation)