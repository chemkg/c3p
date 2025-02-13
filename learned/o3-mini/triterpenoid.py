"""
Classifies: CHEBI:36615 triterpenoid
"""
"""
Classifies: Triterpenoid

Definition: Any terpenoid derived from a triterpene. In our heuristic we first
remove likely sugar (glycoside) moieties using a simple SMARTS, then we identify
the largest fused ring system. We require that the overall aglycone has a total carbon
count between 25 and 70. For the fused core we require it to contain at least 3 rings.
If exactly 3 fused rings are present we demand that the core’s carbon fraction (number
of carbons in the fused core divided by the total carbon count of the aglycone) is at least 50%
(and at least 35% if 4 or more rings). Finally, we check that the raw core carbon count
lies between 15 and 40. (These heuristics are not perfect but are designed to improve
the previous version.)
"""
from rdkit import Chem

def is_triterpenoid(smiles: str):
    """
    Determines if a molecule is a triterpenoid using heuristic rules.
    
    Steps:
      1. Parse the SMILES and count total carbons.
         Reject if the overall carbon count (on the aglycone) is outside a chosen range.
      2. Attempt to remove sugar moieties using a simple pyranose SMARTS.
      3. Identify all rings and then build a graph connecting rings that share atoms.
      4. Determine the largest fused ring cluster.
      5. Count the number of rings in the cluster and gather all atoms in it.
      6. Count the number of carbons in the fused core and compute the ratio relative to the total.
      7. Accept the aglycone if:
           - Total carbons is between 25 and 70,
           - The fused ring system has ≥3 rings,
           - The core carbon count is between 15 and 40, and
           - The ratio of core carbons to total carbons is at least 0.50 when there are only 3 rings,
             or at least 0.35 when there are 4 or more rings.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if classified as a triterpenoid, False otherwise.
        str: Explanation of the decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Count total carbons in the molecule.
    total_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if total_carbons < 25:
        return False, f"Total carbon count ({total_carbons}) is too low to be a triterpenoid"
    if total_carbons > 70:
        return False, f"Total carbon count ({total_carbons}) is too high to be a triterpenoid"
    
    # --- Pre-process: remove likely sugar moieties ---
    # Here we use a simple SMARTS that should catch many common pyranose-like sugar rings.
    sugar_smarts = "[OX2H][CX4]([OX2H])[CX4]([OX2H])[CX4]([OX2H])C"  # a very rough pattern
    sugar_query = Chem.MolFromSmarts(sugar_smarts)
    if sugar_query is not None:
        # Delete substructures that match the sugar pattern.
        aglycone = Chem.DeleteSubstructs(mol, sugar_query)
        aglycone = Chem.RemoveHs(aglycone)
    else:
        aglycone = mol

    # Use the aglycone for further analysis and recalc total carbons using aglycone.
    total_carbons = sum(1 for atom in aglycone.GetAtoms() if atom.GetAtomicNum() == 6)
    if total_carbons < 25 or total_carbons > 70:
        return False, f"After sugar removal, total carbon count ({total_carbons}) is not within acceptable range (25-70)"
    
    # --- Identify fused ring systems ---
    ring_info = aglycone.GetRingInfo()
    rings = ring_info.AtomRings()
    if not rings:
        return False, "No rings detected in the (aglycone) molecule"
    
    # Build a graph where each node is a ring (indexed by its position in 'rings')
    # Two rings are connected if they share at least one atom.
    ring_graph = {i: set() for i in range(len(rings))}
    for i in range(len(rings)):
        for j in range(i+1, len(rings)):
            if set(rings[i]).intersection(rings[j]):
                ring_graph[i].add(j)
                ring_graph[j].add(i)
    
    # Identify connected components (each represents a fused ring system)
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
    
    # Identify the largest fused ring cluster
    largest_component = max(components, key=lambda comp: len(comp))
    num_rings = len(largest_component)
    if num_rings < 3:
        return False, f"The largest fused ring system consists of only {num_rings} ring(s); expected at least 3 for a triterpenoid core"
    
    # Gather all atom indices in the fused ring core.
    core_atom_indices = set()
    for ring_idx in largest_component:
        core_atom_indices.update(rings[ring_idx])
    
    # Count number of carbons in the fused core.
    core_carbons = sum(1 for idx in core_atom_indices if aglycone.GetAtomWithIdx(idx).GetAtomicNum() == 6)
    
    # Check that the core carbon count is in a plausible range.
    if core_carbons < 15 or core_carbons > 40:
        return False, f"Fused ring system core has {core_carbons} carbons, which is not within the acceptable range (15-40)"
    
    # Compute the fraction of (aglycone) carbons that are in the fused core.
    ratio = core_carbons / total_carbons
    # Apply a tighter ratio requirement if the fused system has only 3 rings.
    if num_rings == 3 and ratio < 0.50:
        return False, f"Core carbon fraction ({ratio:.2f}) is too low (need at least 0.50 for a 3-ring system) relative to total carbons ({total_carbons})"
    elif num_rings >= 4 and ratio < 0.35:
        return False, f"Core carbon fraction ({ratio:.2f}) is too low (need at least 0.35 for a {num_rings}-ring system) relative to total carbons ({total_carbons})"
    
    return True, (f"Found triterpene core with {core_carbons} carbons in a fused ring system of {num_rings} rings "
                  f"(total carbons = {total_carbons}, core ratio = {ratio:.2f})")

# Example usage:
if __name__ == "__main__":
    # Test with one of the provided triterpenoid examples (e.g. Neoruscogenin)
    test_smiles = "O1[C@@]2([C@@]([C@@]3([C@]([C@]4([C@](CC3)([C@@]5(C(=CC4)C[C@@H](O)C[C@H]5O)C)[H])[H])(C2)[H])C)([C@@H]([C@]16OCC(CC6)=C)C)[H])[H]"
    result, reason = is_triterpenoid(test_smiles)
    print(result, reason)