"""
Classifies: CHEBI:61655 steroid saponin
"""
#!/usr/bin/env python
"""
Classifies: Steroid saponin
Definition: Any saponin derived from a hydroxysteroid.
That is, the molecule must contain 
  (1) a mostly-carbon steroid tetracyclic core (three 6-membered rings and one 5-membered ring, with at least 17 carbons in the core and at least one exocyclic –OH)
  and 
  (2) at least one sugar ring (a 5- or 6-membered ring containing exactly one ring oxygen and at least one exocyclic –OH)
  that is directly attached to the steroid core.
"""

from rdkit import Chem
from itertools import combinations

def is_steroid_saponin(smiles: str):
    """
    Determines if a molecule is a steroid saponin based on its SMILES string.
    
    The heuristic is as follows:
    1. Parse the SMILES string.
    2. Identify candidate rings (only 5- or 6-membered) that are “mostly carbon”:
         – In a 5-membered ring, at least 4 atoms must be carbon.
         – In a 6-membered ring, at least 5 atoms must be carbon.
    3. Build a graph of rings that are fused (sharing at least 2 atoms) and find connected components.
    4. In each connected component, search for any combination of 4 rings that is composed of one 5-membered ring and three 6-membered rings.
         • The union of atoms in these rings must include at least 17 carbon atoms.
         • Also, at least one carbon in the core must have an exocyclic hydroxyl group attached – meaning it’s connected by a single bond to an oxygen that itself has at least one hydrogen.
    5. Next, search for at least one sugar ring:
         • A sugar ring is defined as a 5- or 6-membered ring containing exactly one ring oxygen,
         • It has at least one exocyclic hydroxyl group,
         • And at least one of its atoms is directly bonded to an atom in the steroid core.
    
    Args:
       smiles (str): A SMILES string of the molecule.
       
    Returns:
       (bool, str): A tuple where the first element is True if the molecule qualifies as a steroid saponin
                    and False otherwise, the second element is a message giving the reasoning.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # For correct detection of exocyclic hydroxyls, add explicit hydrogens.
    molH = Chem.AddHs(mol)
    
    ring_info = mol.GetRingInfo()
    all_rings = ring_info.AtomRings()
    
    # 1. Collect candidate rings: only 5- or 6-membered rings that are mostly carbon.
    candidate_rings = []
    for ring in all_rings:
        if len(ring) not in (5, 6):
            continue
        # Count carbon atoms in the ring.
        nC = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
        if len(ring) == 5 and nC < 4:
            continue
        if len(ring) == 6 and nC < 5:
            continue
        candidate_rings.append(ring)
    
    # 2. Build a graph of fused candidate rings – rings are fused if they share at least 2 atoms.
    ring_graph = {i: set() for i in range(len(candidate_rings))}
    for i in range(len(candidate_rings)):
        for j in range(i+1, len(candidate_rings)):
            if len(set(candidate_rings[i]).intersection(candidate_rings[j])) >= 2:
                ring_graph[i].add(j)
                ring_graph[j].add(i)
    
    # Find connected components in the ring fusion graph using DFS.
    def dfs(node, visited, comp):
        visited.add(node)
        comp.add(node)
        for nei in ring_graph[node]:
            if nei not in visited:
                dfs(nei, visited, comp)
    
    components = []
    visited_nodes = set()
    for node in ring_graph:
        if node not in visited_nodes:
            comp = set()
            dfs(node, visited_nodes, comp)
            components.append(comp)
    
    steroid_core_found = False
    core_atoms = set()
    core_reason = ""
    
    # 3. Search for a candidate steroid core in the connected components.
    # We look for a set of 4 fused rings: one 5-membered and three 6-membered rings.
    for comp in components:
        if len(comp) < 4:
            continue
        # Check every combination of 4 rings.
        for subset in combinations(comp, 4):
            ring_sizes = [len(candidate_rings[i]) for i in subset]
            if sorted(ring_sizes) != [5, 6, 6, 6]:
                continue
            # Combine atoms from these rings; these form the steroid core candidate.
            subset_atoms = set()
            for i in subset:
                subset_atoms.update(candidate_rings[i])
            nC_in_core = sum(1 for idx in subset_atoms if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
            if nC_in_core < 17:
                continue
            # Check for at least one exocyclic hydroxyl group: 
            # look for a neighbor oxygen (bonded by a single bond) attached to a core carbon that is not in the core.
            hydroxyl_found = False
            for idx in subset_atoms:
                atom = mol.GetAtomWithIdx(idx)
                if atom.GetAtomicNum() != 6:
                    continue
                for nbr in atom.GetNeighbors():
                    if nbr.GetIdx() in subset_atoms:
                        continue
                    if nbr.GetAtomicNum() == 8:
                        bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
                        if bond is not None and bond.GetBondType().name == "SINGLE":
                            # Use molH (the molecule with explicit hydrogens) to ensure the oxygen has at least one hydrogen.
                            nbrH = molH.GetAtomWithIdx(nbr.GetIdx())
                            numHs = sum(1 for n in nbrH.GetNeighbors() if n.GetAtomicNum() == 1)
                            if numHs >= 1:
                                hydroxyl_found = True
                                break
                if hydroxyl_found:
                    break
            if hydroxyl_found:
                steroid_core_found = True
                core_atoms = subset_atoms
                core_reason = ("Found fused steroid core (rings: sizes {} with {} core carbons) "
                               "and at least one exocyclic -OH".format(sorted(ring_sizes), nC_in_core))
                break
        if steroid_core_found:
            break
    
    if not steroid_core_found:
        return False, "No steroid (hydroxysteroid) tetracyclic core found"
    
    # 4. Look for a sugar (glycoside) ring attached directly to the steroid core.
    # A sugar ring is a 5- or 6-membered ring that:
    #  - Contains exactly one oxygen atom in the ring,
    #  - Has at least one exocyclic hydroxyl group,
    #  - And at least one of its atoms is directly bonded to an atom in the steroid core.
    sugar_found = False
    sugar_reason = ""
    for ring in all_rings:
        if len(ring) not in (5, 6):
            continue
        oxy_count = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8)
        if oxy_count != 1:
            continue
        # Count exocyclic OH groups attached to ring atoms.
        oh_count = 0
        for idx in ring:
            atom = molH.GetAtomWithIdx(idx)
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring:
                    continue
                if nbr.GetAtomicNum() == 8:
                    bond = molH.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
                    if bond is not None and bond.GetBondType().name == "SINGLE":
                        numHs = sum(1 for n in nbr.GetNeighbors() if n.GetAtomicNum() == 1)
                        if numHs >= 1:
                            oh_count += 1
        if oh_count < 1:
            continue
        # Check if at least one atom in the sugar ring is bonded to an atom in the steroid core.
        attached = False
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in core_atoms:
                    attached = True
                    break
            if attached:
                break
        if attached:
            sugar_found = True
            sugar_reason = ("Sugar ring found with ring oxygen count = 1, exocyclic -OH count = {}, "
                            "attached to steroid core.".format(oh_count))
            break

    if not sugar_found:
        return False, "No glycoside (sugar) ring attached directly to the steroid core found"

    return True, ("Molecule classified as steroid saponin: steroid core detected ({}); sugar moiety detected ({})."
                  .format(core_reason, sugar_reason))

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