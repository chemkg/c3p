"""
Classifies: CHEBI:51963 hopanoid
"""
"""
Classifies: A hopanoid, defined as a triterpenoid based on a hopane skeleton.
A hopanoid is (roughly) expected to have a saturated, fused polycyclic core 
with ≥5 rings (with at least one five‐membered ring and three six‐membered rings)
and a total carbon count in the molecule (or core) near 30. Additional criteria 
are used to reduce false positives (for example, many aromatic rings are not typical).
Note: If the classification is too uncertain, the function may return (None, None).
"""
from rdkit import Chem

def is_hopanoid(smiles: str):
    """
    Determines if a molecule belongs to the hopanoid class based on its SMILES string.
    
    The molecule must have a fused polycyclic core (a connected group of rings, sharing at least
    two atoms per pair) with at least 5 rings, at least one of which is 5-membered and at least three
    that are 6-membered, and the overall molecule should be roughly in the triterpenoid range
    (total carbon count between 25 and 40) with a high sp3 fraction (most carbons are saturated).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a hopanoid, False otherwise.
        str: Explanation for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get overall carbon count in the molecule.
    all_atoms = list(mol.GetAtoms())
    total_carbon = sum(1 for atom in all_atoms if atom.GetAtomicNum() == 6)
    if total_carbon < 25 or total_carbon > 40:
        return False, f"Total carbon count ({total_carbon}) is outside expected triterpenoid range (25-40)"
    
    # Calculate sp3 carbon fraction (hopanoids are highly saturated).
    carbons = [atom for atom in all_atoms if atom.GetAtomicNum() == 6]
    nsp3 = sum(1 for atom in carbons if atom.GetHybridization() == Chem.rdchem.HybridizationType.SP3)
    sp3_fraction = nsp3/len(carbons) if carbons else 0
    if sp3_fraction < 0.7:
        return False, f"Low sp3 fraction ({sp3_fraction:.2f}) among carbons; hopanoids are expected to be highly saturated"
    
    # Get ring information from RDKit. (Each ring is returned as a tuple of atom indices.)
    ring_info = mol.GetRingInfo().AtomRings()
    num_rings = len(ring_info)
    if num_rings < 5:
        return False, f"Found only {num_rings} rings; hopanoid core requires at least 5 fused rings"
    
    # Create sets so we can test “fused” rings (we define rings as fused if they share ≥2 atoms)
    ring_sets = [set(r) for r in ring_info]
    n = len(ring_info)
    adj = {i: set() for i in range(n)}
    for i in range(n):
        for j in range(i+1, n):
            if len(ring_sets[i].intersection(ring_sets[j])) >= 2:
                adj[i].add(j)
                adj[j].add(i)
                
    # Find connected components (clusters of fused rings) via DFS.
    visited = set()
    components = []
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
        components.append(comp)
    
    # Use the largest fused ring cluster for further criteria.
    max_component = max(components, key=lambda comp: len(comp))
    if len(max_component) < 5:
        return False, f"Largest group of fused rings has {len(max_component)} rings; hopanoid core requires at least 5 fused rings"
    
    # Check the ring sizes in the cluster.
    cluster_ring_sizes = [len(ring_info[i]) for i in max_component]
    num_5 = sum(1 for size in cluster_ring_sizes if size == 5)
    num_6 = sum(1 for size in cluster_ring_sizes if size == 6)
    if num_5 < 1:
        return False, "Fused ring cluster does not contain any 5-membered ring; hopanoid core requires at least one"
    if num_6 < 3:
        return False, f"Fused ring cluster contains only {num_6} six-membered rings; at least 3 are required"
    
    # Get all atoms that participate in the fused cluster (union over constituent rings).
    fused_atoms = set()
    for i in max_component:
        fused_atoms |= ring_sets[i]
    
    # Ensure that none of the atoms in the core are aromatic.
    for idx in fused_atoms:
        if mol.GetAtomWithIdx(idx).GetIsAromatic():
            return False, "Fused ring core contains aromatic atoms; hopanoid cores are typically fully alicyclic (non‐aromatic)"
    
    # (Optional) One could try to “recover” additional carbons that belong to the core by expanding
    # one bond out from the fused rings. However, given inconsistent ring perception for hopane-like
    # backbones by RDKit, we use here the union of the rings from the fused cluster.
    fused_core_carbons = sum(1 for idx in fused_atoms if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
    # Note: In some hopanoids the union appears lower than the “true” 30. If fused_core_carbons is very low,
    # we are reluctant to classify as hopanoid.
    if fused_core_carbons < 20:
        return False, f"Fused core appears to contain only {fused_core_carbons} carbon atoms; too few for a hopane skeleton"
    
    # If you reached here, we consider the molecule to meet many of the expected criteria.
    return True, (f"Contains a fused polycyclic structure (largest cluster with {len(max_component)} rings, ring sizes: {cluster_ring_sizes}) "
                  f"and overall molecule has {total_carbon} carbons with high sp3 fraction ({sp3_fraction:.2f}), "
                  f"consistent with a hopanoid skeleton")

# Example usage (for testing):
if __name__ == "__main__":
    # Example: (32R,33R,34R)-bacteriohopanetetrol (a known hopanoid)
    test_smiles = "O[C@@H]([C@H](O)CO)[C@H](O)CC[C@H]([C@@H]1[C@@H]2[C@]([C@@H]3[C@@]([C@]4([C@@H]([C@@]5([C@H](C(CCC5)(C)C)CC4)C)CC3)C)(C)CC2)(C)CC1)C"
    result, reason = is_hopanoid(test_smiles)
    print("Is hopanoid?", result)
    print("Reason:", reason)