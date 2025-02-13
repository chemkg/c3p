"""
Classifies: CHEBI:36615 triterpenoid
"""
"""
Classifies: Triterpenoid
Definition: Any terpenoid derived from a triterpene. The term includes compounds in which the C30 skeleton of the parent triterpene has been rearranged or modified by the removal of one or more skeletal atoms (generally methyl groups).

This heuristic implementation looks for a fused ring‐system of 4 or more rings whose union contains roughly 25–35 carbon atoms (one expects a triterpene core in this range) and that the overall molecule contains at least 27 carbons. Note that many triterpenoids are decorated with sugars or extra oxygen functional groups so our test is only approximate.
"""

from rdkit import Chem

def is_triterpenoid(smiles: str):
    """
    Determines if the molecule is a triterpenoid using heuristic rules.
    
    Heuristics used:
      1. The SMILES string must parse into a molecule.
      2. The molecule must contain rings.
      3. We identify fused ring systems by checking rings (via GetRingInfo)
         that share atoms. Triterpene cores are usually built of 4 or more fused rings.
      4. The fused ring system (largest connected ring component) is assumed to be the triterpene
         core. We count the number of carbon atoms in this core – typically between 25 and 35.
      5. The overall molecule should have at least 27 carbon atoms.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a triterpenoid, False otherwise.
        str: Reason for the classification decision.
    """
    # Parse SMILES into molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Calculate total carbon count in molecule
    total_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if total_carbons < 27:
        return False, f"Total carbon count ({total_carbons}) is too low to be a triterpenoid"
    
    # Retrieve ring information
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    if not rings:
        return False, "No rings detected in molecule"
    
    # Build a connectivity (graph) between rings: each node is a ring (indexed by its order in 'rings')
    # Two rings are connected if they share at least one common atom.
    ring_graph = {i: set() for i in range(len(rings))}
    for i in range(len(rings)):
        for j in range(i+1, len(rings)):
            if set(rings[i]).intersection(rings[j]):
                ring_graph[i].add(j)
                ring_graph[j].add(i)
    
    # Find all connected components in the ring graph using a simple DFS
    visited = set()
    components = []
    
    for i in range(len(rings)):
        if i not in visited:
            stack = [i]
            comp = set()
            while stack:
                node = stack.pop()
                if node not in visited:
                    visited.add(node)
                    comp.add(node)
                    stack.extend(ring_graph[node] - visited)
            components.append(comp)
    
    # Identify the largest connected ring cluster (core)
    largest_component = max(components, key=lambda comp: len(comp))
    # Get combined set of atoms in the largest fused ring system 
    core_atom_indices = set()
    for ring_idx in largest_component:
        core_atom_indices.update(rings[ring_idx])
    # Count carbon atoms in the core
    core_carbons = sum(1 for idx in core_atom_indices if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
    
    # Heuristic: a triterpene core normally has about 30 carbons,
    # but rearrangements (loss or migration of methyl groups) might lower that number.
    if core_carbons < 25 or core_carbons > 35:
        return False, f"Fused ring system core has {core_carbons} carbons, which is not within the expected range (25-35) for a triterpene core"
    
    # Also require that the fused ring system is built from at least 4 rings.
    if len(largest_component) < 4:
        return False, f"The largest fused ring system consists of only {len(largest_component)} rings; expected at least 4 for a triterpene core"
    
    return True, f"Found triterpene core with {core_carbons} carbons in a fused ring system of {len(largest_component)} rings (total carbons = {total_carbons})"
    
# Example usage (if needed for testing):
if __name__ == "__main__":
    # Example triterpenoid: Neoruscogenin (from the prompt)
    test_smiles = "O1[C@@]2([C@@]([C@@]3([C@]([C@]4([C@](CC3)([C@@]5(C(=CC4)C[C@@H](O)C[C@H]5O)C)[H])[H])(C2)[H])C)([C@@H]([C@]16OCC(CC6)=C)C)[H])[H]"
    result, reason = is_triterpenoid(test_smiles)
    print(result, reason)