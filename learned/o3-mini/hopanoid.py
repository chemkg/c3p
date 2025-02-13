"""
Classifies: CHEBI:51963 hopanoid
"""
"""
Classifies: A hopanoid, defined as a triterpenoid based on a hopane skeleton.
We use the heuristic that a hopanoid should have a fused polycyclic core where at
least 5 rings are interconnected (sharing at least 2 atoms per adjacent ring) and
a sufficient number of carbon atoms (usually around 30 for a triterpenoid).
"""
from rdkit import Chem

def is_hopanoid(smiles: str):
    """
    Determines if a molecule belongs to the hopanoid class based on its SMILES string.
    A hopanoid is defined as a triterpenoid based on a hopane skeleton – in other words,
    the molecule must have a fused polycyclic core with (at least) five rings that share 
    bonds and a high carbon content.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is recognized as a hopanoid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES into RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Obtain ring information from the molecule.
    ring_info = mol.GetRingInfo().AtomRings()
    num_rings = len(ring_info)
    if num_rings < 5:
        return False, f"Found only {num_rings} rings; hopanoid core requires at least 5 fused rings"
    
    # Build a graph among rings (nodes = rings); two rings are considered fused if they share
    # two or more atoms.
    n = len(ring_info)
    # Create a list of sets for each ring's atom indices for easy intersection tests.
    ring_sets = [set(ring) for ring in ring_info]
    # Build adjacency list: for each ring, list indices of rings fused to it.
    adj = {i: set() for i in range(n)}
    for i in range(n):
        for j in range(i+1, n):
            if len(ring_sets[i].intersection(ring_sets[j])) >= 2:
                adj[i].add(j)
                adj[j].add(i)

    # Now find connected components (clusters of fused rings) in the graph.
    visited = set()
    components = []
    for i in range(n):
        if i in visited:
            continue
        # perform DFS for component starting at ring i
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
    
    # Check if any connected component has at least 5 rings.
    max_component_size = max(len(comp) for comp in components)
    if max_component_size < 5:
        return False, (f"Largest group of fused rings has {max_component_size} rings; "
                       f"hopanoid core requires at least 5 fused rings")
    
    # Check the carbon count as hopanoids (triterpenoids) are mainly built from carbons.
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 27:
        return False, f"Found only {carbon_count} carbon atoms; typical hopanoid triterpenoids have around 30 carbons"
    
    return True, f"Contains fused polycyclic structure (largest cluster size: {max_component_size} rings) and {carbon_count} carbons, consistent with a hopanoid skeleton"

# Example usage (for testing):
if __name__ == "__main__":
    # Example: (32R,33R,34S)-bacteriohopanetetrol
    test_smiles = "O[C@@H]([C@@H](O)CO)[C@H](O)CC[C@H]([C@@H]1[C@@H]2[C@]([C@@H]3[C@@]([C@]4([C@@H]([C@@]5([C@H](C(CCC5)(C)C)CC4)C)CC3)C)(C)CC2)(C)CC1)C"
    result, reason = is_hopanoid(test_smiles)
    print("Is hopanoid?", result)
    print("Reason:", reason)