"""
Classifies: CHEBI:51963 hopanoid
"""
"""
Classifies: A hopanoid, defined as a triterpenoid based on a hopane skeleton.
A true hopanoid should have a fused (≥5 rings) polycyclic core with a characteristic 
6–6–6–6–5 ring arrangement (i.e. at least one 5-membered ring and at least 3 six‐membered rings 
in the fused system) and a sufficient number of carbon atoms in the core (around 30).
"""
from rdkit import Chem

def is_hopanoid(smiles: str):
    """
    Determines if a molecule belongs to the hopanoid class based on its SMILES string.
    
    A hopanoid is defined as a triterpenoid based on a hopane skeleton – that is, the molecule
    should possess a fused polycyclic core of at least 5 rings. In hopanoids the pentacyclic skeleton 
    tends to have a specific 6–6–6–6–5 ring arrangement (i.e. one of the rings is 5-membered, while 
    at least three rings are 6-membered) and the fused core should contain roughly 30 carbon atoms.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is recognized as a hopanoid, False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get ring information: each ring is a tuple of atom indices.
    ring_info = mol.GetRingInfo().AtomRings()
    num_rings = len(ring_info)
    
    if num_rings < 5:
        return False, f"Found only {num_rings} rings; hopanoid core requires at least 5 fused rings"
    
    # For each ring, get a set of atom indices.
    ring_sets = [set(ring) for ring in ring_info]
    
    # Build an adjacency graph between rings: two rings are fused if they share >=2 atoms.
    n = len(ring_info)
    adj = {i: set() for i in range(n)}
    for i in range(n):
        for j in range(i+1, n):
            if len(ring_sets[i].intersection(ring_sets[j])) >= 2:
                adj[i].add(j)
                adj[j].add(i)
    
    # Find connected components (clusters of fused rings) using DFS.
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
    
    # Identify the largest fused cluster.
    max_component = max(components, key=lambda comp: len(comp))
    if len(max_component) < 5:
        return False, (f"Largest group of fused rings has {len(max_component)} rings; "
                       f"hopanoid core requires at least 5 fused rings")
    
    # Now, for the rings in the largest fused cluster, check ring sizes.
    # We extract the list of ring sizes.
    ring_sizes = [len(ring_info[i]) for i in max_component]
    
    # Check for presence of at least one 5-membered ring and at least three 6-membered rings.
    num_5 = sum(1 for size in ring_sizes if size == 5)
    num_6 = sum(1 for size in ring_sizes if size == 6)
    if num_5 < 1:
        return False, f"Largest fused ring cluster does not contain any 5-membered ring; hopanoid core requires a 5-membered ring"
    if num_6 < 3:
        return False, f"Largest fused ring cluster contains only {num_6} six-membered rings; hopanoid core requires at least 3 six-membered rings"
    
    # Next, measure the carbon count in the fused core: take the union of atoms in the largest cluster.
    fused_atom_indices = set()
    for i in max_component:
        fused_atom_indices |= ring_sets[i]
    
    # Count only carbon atoms in the fused core.
    carbon_count = sum(1 for idx in fused_atom_indices if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
    if carbon_count < 30:
        return False, f"Fused core contains only {carbon_count} carbon atoms; expected at least 30 carbons in a hopanoid skeleton"
    
    return True, (f"Contains fused polycyclic structure (largest cluster of {len(max_component)} rings with "
                  f"ring sizes {ring_sizes}) and fused core has {carbon_count} carbons, consistent with a hopanoid skeleton")

# Example usage (for testing):
if __name__ == "__main__":
    # Example: (32R,33R,34R)-bacteriohopanetetrol (a known hopanoid)
    test_smiles = "O[C@@H]([C@H](O)CO)[C@H](O)CC[C@H]([C@@H]1[C@@H]2[C@]([C@@H]3[C@@]([C@]4([C@@H]([C@@]5([C@H](C(CCC5)(C)C)CC4)C)CC3)C)(C)CC2)(C)CC1)C"
    result, reason = is_hopanoid(test_smiles)
    print("Is hopanoid?", result)
    print("Reason:", reason)