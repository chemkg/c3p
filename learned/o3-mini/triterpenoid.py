"""
Classifies: CHEBI:36615 triterpenoid
"""
"""
Classifies: Triterpenoid

Definition: Any terpenoid derived from a triterpene. In our heuristic we first
remove likely sugar (glycoside) moieties using a simple SMARTS. Then, after verifying
that the (aglycone) total carbon count is between 25 and 70, we identify the largest fused 
ring system. We require that the fused core consists of at least 3 rings, that the core's 
carbon count is between 15 and 40, and that the ratio of core carbons to aglycone carbons
is at least 0.50 (if 3 rings) or 0.35 (if 4 or more rings). Additionally, we now require 
that a sufficient fraction of core carbons are sp3-hybridized (≥50%), as triterpenoids 
usually have saturated (or partially saturated) fused cores.
"""
from rdkit import Chem

def is_triterpenoid(smiles: str):
    """
    Determines if a molecule is a triterpenoid using heuristic rules.
    
    Steps:
      1. Parse the SMILES and count total carbons.
         Reject if the overall carbon count (on the aglycone) is outside 25--70.
      2. Remove likely sugar moieties using a simple pyranose-like SMARTS.
      3. Recalculate aglycone total carbons.
      4. Identify rings and build a graph connecting rings that share atoms.
      5. Find the largest connected (fused) ring system.
      6. Count carbons in this fused core.
      7. Compute the ratio of core carbons to total aglycone carbons. For a 3-ring
         system the ratio should be at least 0.50; for systems with 4 or more rings, at least 0.35.
      8. Additionally, count how many core carbons are sp3-hybridized. In a typical
         triterpenoid, at least 50% of core carbons should be sp3.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if classified as a triterpenoid, False otherwise.
        str: Explanation of the classification decision.
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
    # We use a simple (rough) SMARTS pattern to match common pyranose-like sugars.
    sugar_smarts = "[OX2H][CX4]([OX2H])[CX4]([OX2H])[CX4]([OX2H])C"
    sugar_query = Chem.MolFromSmarts(sugar_smarts)
    if sugar_query is not None:
        aglycone = Chem.DeleteSubstructs(mol, sugar_query)
        aglycone = Chem.RemoveHs(aglycone)
    else:
        aglycone = mol

    # Recalculate total carbons on the aglycone.
    total_carbons = sum(1 for atom in aglycone.GetAtoms() if atom.GetAtomicNum() == 6)
    if total_carbons < 25 or total_carbons > 70:
        return False, f"After sugar removal, total carbon count ({total_carbons}) is not within acceptable range (25-70)"

    # --- Identify fused ring systems ---
    ring_info = aglycone.GetRingInfo()
    rings = ring_info.AtomRings()
    if not rings:
        return False, "No rings detected in the (aglycone) molecule"

    # Build a graph: nodes are rings (indexed by position in 'rings'), connected if they share at least one atom.
    ring_graph = {i: set() for i in range(len(rings))}
    for i in range(len(rings)):
        for j in range(i+1, len(rings)):
            if set(rings[i]).intersection(rings[j]):
                ring_graph[i].add(j)
                ring_graph[j].add(i)

    # Identify connected components (each component is a fused ring system)
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

    # Identify the largest fused ring cluster.
    largest_component = max(components, key=lambda comp: len(comp))
    num_rings = len(largest_component)
    if num_rings < 3:
        return False, f"The largest fused ring system consists of only {num_rings} ring(s); expected at least 3 for a triterpenoid core"

    # Gather all unique atom indices in the fused ring core.
    core_atom_indices = set()
    for ring_idx in largest_component:
        core_atom_indices.update(rings[ring_idx])
    
    # Count number of carbons in the fused core.
    core_carbons = sum(1 for idx in core_atom_indices if aglycone.GetAtomWithIdx(idx).GetAtomicNum() == 6)
    if core_carbons < 15 or core_carbons > 40:
        return False, f"Fused ring system core has {core_carbons} carbons, which is not within the acceptable range (15-40)"
    
    # Compute the ratio of core carbons to total aglycone carbons.
    core_ratio = core_carbons / total_carbons
    if num_rings == 3 and core_ratio < 0.50:
        return False, f"Core carbon fraction ({core_ratio:.2f}) is too low (need at least 0.50 for a 3-ring system) relative to total carbons ({total_carbons})"
    elif num_rings >= 4 and core_ratio < 0.35:
        return False, f"Core carbon fraction ({core_ratio:.2f}) is too low (need at least 0.35 for a {num_rings}-ring system) relative to total carbons ({total_carbons})"
    
    # --- Additional Saturation Check: Ensure the fused core is sufficiently saturated ---
    # Count how many of the core carbons are sp3 hybridized.
    sp3_core = 0
    for idx in core_atom_indices:
        atom = aglycone.GetAtomWithIdx(idx)
        if atom.GetAtomicNum() == 6 and atom.GetHybridization() == Chem.rdchem.HybridizationType.SP3:
            sp3_core += 1
    sp3_fraction = sp3_core / core_carbons
    if sp3_fraction < 0.50:
        return False, (f"Fused ring core is too unsaturated: only {sp3_fraction:.2f} of its {core_carbons} carbons are sp3; "
                       "expected ≥0.50 for a triterpenoid")
    
    return True, (f"Found triterpene core with {core_carbons} carbons in a fused ring system of {num_rings} rings "
                  f"(total aglycone carbons = {total_carbons}, core ratio = {core_ratio:.2f}, sp3 fraction = {sp3_fraction:.2f})")

# Example usage:
if __name__ == "__main__":
    # Test one example: Neoruscogenin
    test_smiles = "O1[C@@]2([C@@]([C@@]3([C@]([C@]4([C@](CC3)([C@@]5(C(=CC4)C[C@@H](O)C[C@H]5O)C)[H])[H])(C2)[H])C)([C@@H]([C@]16OCC(CC6)=C)C)[H])[H]"
    result, reason = is_triterpenoid(test_smiles)
    print(result, reason)