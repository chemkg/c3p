"""
Classifies: CHEBI:61655 steroid saponin
"""
#!/usr/bin/env python
"""
Classifies: Steroid saponin
Definition: Any saponin derived from a hydroxysteroid.
That is, a molecule must contain (1) a steroid (hydroxysteroid) tetracyclic core and 
(2) at least one sugar moiety (a glycoside ring) attached as a substituent.
Our approach uses ring analysis rather than a very specific SMARTS pattern.
"""

from rdkit import Chem
from rdkit.Chem import rdmolops

def is_steroid_saponin(smiles: str):
    """
    Determines if a molecule is a steroid saponin based on its SMILES string.
    
    Our strategy (heuristic):
      1. Parse the SMILES string.
      2. Identify candidate rings that could form a steroid core.
         We look for four fused rings that are made entirely of carbons.
         Furthermore, we require that within that fused set there are three 6-membered rings
         and one 5-membered ring.
      3. Verify that at least one carbon in these rings has an exocyclic –OH group (to meet the
         “hydroxysteroid” requirement).
      4. Look for at least one sugar ring – here defined as a 5- or 6-membered ring that 
         contains exactly one oxygen atom in the ring and at least two exocyclic hydroxyl groups.
    
    Args:
      smiles (str): SMILES string of the molecule.
      
    Returns:
      (bool, str): A tuple where the first element is True if the molecule is classified as a 
                   steroid saponin, and False otherwise; the second element is a message explaining 
                   the reasoning.
    """
    # Parse the SMILES; if invalid, return immediately.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # For sugar detection we may need explicit hydrogens.
    mol_with_H = Chem.AddHs(mol)
    
    # Get ring information from the molecule.
    ring_info = mol.GetRingInfo()
    all_rings = ring_info.AtomRings()
    
    # Step 1 — Look for candidate steroid rings:
    # We assume that steroid rings are composed entirely of carbons.
    candidate_rings = []
    for ring in all_rings:
        # Only consider rings that are either 5 or 6 members.
        if len(ring) not in (5, 6):
            continue
        # Check that every atom in the ring is carbon.
        if all(mol.GetAtomWithIdx(idx).GetAtomicNum() == 6 for idx in ring):
            candidate_rings.append(ring)
    
    # Build a graph over candidate rings. Two rings are "fused" if they share at least two atoms.
    ring_graph = {i: set() for i in range(len(candidate_rings))}
    for i in range(len(candidate_rings)):
        for j in range(i+1, len(candidate_rings)):
            shared = set(candidate_rings[i]).intersection(set(candidate_rings[j]))
            if len(shared) >= 2:
                ring_graph[i].add(j)
                ring_graph[j].add(i)
    
    # Now find connected components (clusters) from this graph.
    def dfs(node, visited, comp):
        visited.add(node)
        comp.add(node)
        for neighbor in ring_graph[node]:
            if neighbor not in visited:
                dfs(neighbor, visited, comp)
    
    components = []
    visited = set()
    for node in ring_graph:
        if node not in visited:
            comp = set()
            dfs(node, visited, comp)
            components.append(comp)
    
    steroid_core_found = False
    reason_core = ""
    # Check each component to see if any has exactly 4 rings with the expected size: one 5-membered and three 6-membered.
    for comp in components:
        if len(comp) < 4:
            continue
        # We iterate over all subsets of 4 rings from this component:
        from itertools import combinations
        for subset in combinations(comp, 4):
            ring_sizes = [len(candidate_rings[i]) for i in subset]
            if sorted(ring_sizes) == [5, 6, 6, 6]:
                # Found a candidate steroid nucleus.
                # Now check for at least one hydroxyl (–OH) attached to one of the carbons in these rings.
                hydroxyl_found = False
                # Get all atom indices in these candidate rings.
                core_atoms = set()
                for i in subset:
                    core_atoms.update(candidate_rings[i])
                for idx in core_atoms:
                    atom = mol.GetAtomWithIdx(idx)
                    if atom.GetAtomicNum() != 6:
                        continue
                    # Check neighboring atoms that are not in the core.
                    for nbr in atom.GetNeighbors():
                        if nbr.GetIdx() in core_atoms:
                            continue
                        # Check if the neighbor is oxygen connected via a single bond.
                        if nbr.GetAtomicNum() == 8:
                            bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
                            if bond is not None and bond.GetBondType().name == "SINGLE":
                                hydroxyl_found = True
                                break
                    if hydroxyl_found:
                        break
                if hydroxyl_found:
                    steroid_core_found = True
                    break
        if steroid_core_found:
            break

    if not steroid_core_found:
        return False, "No steroid (hydroxysteroid) tetracyclic core found"
    
    # Step 2 — Look for at least one sugar ring.
    # We use a similar ring search but now we look for rings that are 5 or 6 members,
    # have exactly one ring oxygen (atomic number 8), and at least two exocyclic hydroxyl groups.
    sugar_found = False
    for ring in all_rings:
        if len(ring) not in (5, 6):
            continue
        # Count atoms in the ring that are oxygen. (Note: we use the original mol.)
        oxy_in_ring = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8)
        if oxy_in_ring != 1:
            continue
        # Now count exocyclic –OH groups attached to ring atoms.
        oh_count = 0
        for idx in ring:
            atom = mol_with_H.GetAtomWithIdx(idx)
            # We expect sugar rings to consist mainly of carbon.
            if atom.GetAtomicNum() != 6:
                continue
            for nbr in atom.GetNeighbors():
                # If the neighbor is not in this ring...
                if nbr.GetIdx() in ring:
                    continue
                # Check if neighboring oxygen is part of an –OH.
                if nbr.GetAtomicNum() == 8:
                    bond = mol_with_H.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
                    if bond is not None and bond.GetBondType().name == "SINGLE":
                        # Further check: the oxygen should have at least one hydrogen.
                        # (Counting explicit H atoms attached.)
                        num_H = sum(1 for h in nbr.GetNeighbors() if h.GetAtomicNum() == 1)
                        if num_H >= 1:
                            oh_count += 1
        if oh_count >= 2:
            sugar_found = True
            break

    if not sugar_found:
        return False, "No sugar ring (glycoside) found"
    
    return True, "Molecule contains a hydroxysteroid tetracyclic core with at least one sugar ring attached"

# Example usage (for testing):
if __name__ == "__main__":
    test_smiles = [
        # Jurubine
        "O1[C@@]2([C@@]([C@@]3([C@]([C@]4([C@@]([C@@]5([C@@](CC4)(C[C@@H](N)CC5)[H])C)(CC3)[H])[H])(C2)[H])C)([C@@H]([C@@]1(O)CC[C@H](CO[C@@H]6O[C@@H]([C@@H](O)C(O)C6O)CO)C)C)[H])[H]",
        # Asparagoside B
        "O1[C@@]2([C@@]([C@@]3([C@]([C@]4([C@@]([C@@]5([C@](CC4)(C[C@@H](O)CC5)[H])C)(CC3)[H])[H])(C2)[H])C)([C@@H]([C@@]1(O)CC[C@@H](CO[C@@H]6O[C@@H]([C@@H](O)[C@H](O)[C@H]6O)CO)C)C)[H])[H]",
        # ginsenoside Re
        "C[C@@H]1O[C@@H](O[C@@H]2[C@@H](O)[C@H](O)[C@@H](CO)O[C@H]2O[C@H]2C[C@]3(C)[C@H](C[C@@H](O)[C@@H]4[C@H](CC[C@@]34C)[C@](C)(CCC=C(C)C)O[C@@H]3O[C@H](CO)[C@@H](O)[C@H](O)[C@H]3O)[C@@]3(C)CC[C@H](O)C(C)(C)[C@H]23)[C@H](O)[C@H](O)[C@H]1O",
        # 17-beta-Estradiol glucuronide
        "O([C@@H]1[C@@]2(C(C3C(CC2)C4=C(CC3)C=C(O)C=C4)CC1)C)[C@@H]5O[C@@H]([C@@H](O)[C@H](O)[C@H]5O)C(O)=O"
    ]
    for smi in test_smiles:
        result, reason = is_steroid_saponin(smi)
        print("SMILES:", smi)
        print("Is steroid saponin?", result)
        print("Reason:", reason)
        print("-----")