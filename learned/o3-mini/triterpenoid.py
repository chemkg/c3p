"""
Classifies: CHEBI:36615 triterpenoid
"""
"""
Classifies: Triterpenoid
Definition: Any terpenoid derived from a triterpene. The term includes compounds in which the C30 
skeleton of the parent triterpene has been rearranged or modified by the removal of one or more 
skeletal atoms (generally methyl groups). Because many triterpenoids bear sugars or other decorations, 
this heuristic:
  (1) removes common sugar motifs,
  (2) finds the largest fused ring system, and
  (3) checks that its carbon count and connectivity approximate that expected from a modified 
      triterpene.
Heuristics used:
   - Molecule must parse and contain rings.
   - Optionally remove sugar moieties (detected by a simple pyranose SMARTS).
   - The largest connected (fused) ring cluster must contain at least 4 rings.
   - The “core” carbon count from that cluster is accepted if it lies between a more relaxed range 
     (here we allow between ~15 and 40 carbons, to cover modified/rearranged cases).
   - Additionally, we require that the fused core represents a sufficiently high fraction 
     of the total carbons (here at least 35%).
Note: This is a heuristic and will mis‐classify some genuine triterpenoids and may accept some 
non-triterpenoids.
"""
from rdkit import Chem

def is_triterpenoid(smiles: str):
    """
    Determines if the molecule is a triterpenoid using heuristic rules.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a triterpenoid, False otherwise.
        str: Explanation of the result.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count total carbons in the molecule.
    total_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if total_carbons < 27:
        return False, f"Total carbon count ({total_carbons}) is too low to be a triterpenoid"

    # --- Pre-process: try to remove sugar moieties ---
    # (common sugars often are 6-membered rings with 5 carbons/1 oxygen; this simple SMARTS 
    # may catch many glycosides)
    sugar_smarts = "[C;R6][O;R6][C;R6]([O;R6])[C;R6]([O;R6])[C;R6]1" 
    sugar_mol = Chem.MolFromSmarts(sugar_smarts)
    if sugar_mol is not None:
        aglycone = Chem.DeleteSubstructs(mol, sugar_mol) 
        # Chem.DeleteSubstructs returns a Mol that might be disconnected.
        aglycone = Chem.RemoveHs(aglycone)
    else:
        aglycone = mol

    # Use the aglycone for fused ring system analysis.
    target = aglycone

    # Look for rings in the target molecule.
    ring_info = target.GetRingInfo()
    rings = ring_info.AtomRings()
    if not rings:
        return False, "No rings detected in the (aglycone) molecule"

    # Build a graph where each node is a ring (by its index in rings)
    # Two rings are adjacent if they share at least one atom.
    ring_graph = {i: set() for i in range(len(rings))}
    for i in range(len(rings)):
        for j in range(i+1, len(rings)):
            if set(rings[i]).intersection(rings[j]):
                ring_graph[i].add(j)
                ring_graph[j].add(i)
                
    # Identify connected components (each is a fused ring system)
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
    
    # Identify largest fused ring cluster (by number of rings).
    largest_component = max(components, key=lambda comp: len(comp))
    
    # Require that the fused ring system has at least 4 rings.
    if len(largest_component) < 4:
        return False, f"The largest fused ring system consists of only {len(largest_component)} rings; expected at least 4 for a triterpene core"
    
    # Gather all atom indices in the fused ring core.
    core_atom_indices = set()
    for ring_idx in largest_component:
        core_atom_indices.update(rings[ring_idx])
    
    # Count number of carbons in the fused core.
    core_carbons = sum(1 for idx in core_atom_indices if target.GetAtomWithIdx(idx).GetAtomicNum() == 6)
    core_total_atoms = len(core_atom_indices)
    
    # Heuristic: even when modified, the core (aglycone) should have a carbon count not too far from 30.
    # Here we relax the lower bound to 15 and upper bound to 40.
    if core_carbons < 15 or core_carbons > 40:
        return False, f"Fused ring system core has {core_carbons} carbons, which is not within the acceptable range (15-40) for a triterpene core"
    
    # Also require that the fused core represents a decent fraction of the overall carbons.
    ratio = core_carbons / total_carbons
    if ratio < 0.35:
        return False, f"Core carbon fraction ({ratio:.2f}) is too low relative to total carbons ({total_carbons})"
    
    return True, f"Found triterpene core with {core_carbons} carbons in a fused ring system of {len(largest_component)} rings (total carbons = {total_carbons}, core ratio = {ratio:.2f})"

# Example usage:
if __name__ == "__main__":
    # You can test with one of the provided examples.
    test_smiles = "O1[C@@]2([C@@]([C@@]3([C@]([C@]4([C@](CC3)([C@@]5(C(=CC4)C[C@@H](O)C[C@H]5O)C)[H])[H])(C2)[H])C)([C@@H]([C@]16OCC(CC6)=C)C)[H])[H]"
    result, reason = is_triterpenoid(test_smiles)
    print(result, reason)