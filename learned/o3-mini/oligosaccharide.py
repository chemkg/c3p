"""
Classifies: CHEBI:50699 oligosaccharide
"""
"""
Classifies: Oligosaccharide
Definition: A compound in which monosaccharide units are joined by glycosidic linkages.
Heuristic improvements in this version:
  1. Candidate sugar rings are identified as either 5‐membered (furanose) or 6‐membered (pyranose) rings
     that contain exactly one ring oxygen (with allowance for one nitrogen substituting a carbon) and must be
     fully saturated (all bonds inside the ring are single). Additionally, at least two ring carbons are required
     to have exocyclic hydroxyl groups (i.e. an oxygen neighbor that carries at least one hydrogen).
  2. At least two such sugar rings must be found and they must be connected by a glycosidic linkage.
     Two rings are declared “connected” if any atom in one is directly bond‐linked to an atom in the other, or
     if an atom in one ring is attached to an oxygen (with degree 2) that in turn is bonded to an atom in the other ring.
  3. Finally, the sugar rings should account for at least 50% of the molecule’s heavy atoms.
Note: This is a heuristic approach that will not cover all edge cases.
"""

from rdkit import Chem

def is_oligosaccharide(smiles: str):
    """
    Determines if a molecule is an oligosaccharide based on its SMILES string.
    Uses several heuristics:
      - Identify candidate sugar rings of size 5 or 6 that are saturated (all bonds in the ring are single)
        and contain exactly one oxygen (with an allowance for one nitrogen substitution). 
      - Require that at least two of the ring carbons have an exocyclic hydroxyl group.
      - Require that at least two candidate sugar rings are found and that these rings are connected
        via either a direct bond or a bridging (low-degree) oxygen.
      - Verify that the atoms in the rings make up at least 50% of the heavy atoms.
      
    Args:
       smiles (str): SMILES string of the molecule.
       
    Returns:
       bool: True if the molecule is classified as an oligosaccharide.
       str: Reason for the classification.
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get all heavy atoms (atomic number > 1)
    heavy_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1]
    total_heavy = len(heavy_atoms)
    
    ring_info = mol.GetRingInfo().AtomRings()
    sugar_rings = []  # Each candidate is stored as a set of atom indices
    
    # Loop over each ring in the molecule
    for ring in ring_info:
        if len(ring) not in (5, 6):
            continue
        
        # Check saturation: ensure that every bond connecting two atoms in the ring is a single bond.
        ring_is_saturated = True
        ring_set = set(ring)
        # examine all bonds in the molecule that have both endpoints in the ring
        for bond in mol.GetBonds():
            a_idx = bond.GetBeginAtomIdx()
            b_idx = bond.GetEndAtomIdx()
            if a_idx in ring_set and b_idx in ring_set:
                if bond.GetBondType() != Chem.BondType.SINGLE:
                    ring_is_saturated = False
                    break
        if not ring_is_saturated:
            continue
        
        # Count ring atoms by atomic number: oxygen (8), carbon (6) and nitrogen (7)
        o_count, c_count, n_count = 0, 0, 0
        for idx in ring:
            atomic_num = mol.GetAtomWithIdx(idx).GetAtomicNum()
            if atomic_num == 8:
                o_count += 1
            elif atomic_num == 6:
                c_count += 1
            elif atomic_num == 7:
                n_count += 1
        
        # For a typical sugar ring: exactly one oxygen.
        # For a 6-membered ring (pyranose): expect 1 oxygen and either 5 carbons or 4 carbons plus 1 nitrogen.
        # For a 5-membered ring (furanose): expect 1 oxygen and either 4 carbons or 3 carbons plus 1 nitrogen.
        if len(ring) == 6:
            if not (o_count == 1 and (c_count == 5 or (c_count == 4 and n_count == 1))):
                continue
        elif len(ring) == 5:
            if not (o_count == 1 and (c_count == 4 or (c_count == 3 and n_count == 1))):
                continue
        
        # Additional check: require that at least 2 ring carbons have exocyclic -OH groups.
        hydroxyl_count = 0
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            # Only check for carbons (typical sugar ring atoms except the ring oxygen)
            if atom.GetAtomicNum() != 6:
                continue
            # Look among neighbors that are not in the ring
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring_set:
                    continue
                # Check if neighbor is oxygen and (heuristically) has at least one attached hydrogen.
                if nbr.GetAtomicNum() == 8 and nbr.GetTotalNumHs() > 0:
                    hydroxyl_count += 1
                    break  # count each ring carbon at most once
        if hydroxyl_count < 2:
            continue
        
        # If all checks pass, add this candidate ring.
        sugar_rings.append(ring_set)
    
    if len(sugar_rings) < 2:
        return False, f"Fewer than 2 candidate sugar rings detected (found {len(sugar_rings)})."
    
    # Build connectivity graph among sugar rings.
    # Two rings are "connected" if:
    #  (a) there is a bond directly between an atom in one ring and an atom in the other, OR
    #  (b) an atom in one ring is attached to an oxygen (with degree == 2) that in turn is bonded to an atom in the other ring.
    n = len(sugar_rings)
    graph = {i: set() for i in range(n)}
    for i in range(n):
        for j in range(i + 1, n):
            connected = False
            # Option (a): direct bond between any two ring atoms.
            for bond in mol.GetBonds():
                a_idx = bond.GetBeginAtomIdx()
                b_idx = bond.GetEndAtomIdx()
                if ((a_idx in sugar_rings[i] and b_idx in sugar_rings[j]) or 
                    (a_idx in sugar_rings[j] and b_idx in sugar_rings[i])):
                    connected = True
                    break
            # Option (b): bridging oxygen with degree 2.
            if not connected:
                for a_idx in sugar_rings[i]:
                    atom_a = mol.GetAtomWithIdx(a_idx)
                    for nbr in atom_a.GetNeighbors():
                        # Only count a bridging oxygen if it is not in the same ring, is oxygen, and its degree is 2.
                        if nbr.GetIdx() in sugar_rings[i]:
                            continue
                        if nbr.GetAtomicNum() != 8 or nbr.GetDegree() != 2:
                            continue
                        # Check if the other neighbor of this oxygen is in ring j.
                        for nbr2 in nbr.GetNeighbors():
                            if nbr2.GetIdx() == a_idx:
                                continue
                            if nbr2.GetIdx() in sugar_rings[j]:
                                connected = True
                                break
                        if connected:
                            break
                    if connected:
                        break
            if connected:
                graph[i].add(j)
                graph[j].add(i)
    
    # Check if all candidate sugar rings form a single connected component.
    visited = set()
    stack = [0]
    while stack:
        node = stack.pop()
        if node in visited:
            continue
        visited.add(node)
        for nbr in graph[node]:
            if nbr not in visited:
                stack.append(nbr)
    if len(visited) != n:
        return False, "Candidate sugar rings are detected but they are not all connected via glycosidic bonds."
    
    # Finally, count the heavy atoms that are part of any sugar ring.
    sugar_atoms = set()
    for ring in sugar_rings:
        sugar_atoms.update(ring)
    fraction = len(sugar_atoms) / total_heavy if total_heavy else 0
    if fraction < 0.5:
        return False, f"Sugar rings detected but they only constitute {fraction*100:.1f}% of heavy atoms; likely a glycoside rather than an oligosaccharide."
    
    return True, "Molecule has multiple connected, saturated sugar rings with characteristic hydroxyls and these rings account for the majority of heavy atoms; likely an oligosaccharide."

# Example usage:
if __name__ == '__main__':
    # Example: alpha-mannobiose (a disaccharide)
    test_smiles = "OC[C@H]1O[C@@H](O[C@@H]2[C@@H](CO)O[C@H](O)[C@@H](O)[C@H]2O)[C@@H](O)[C@@H](O)[C@@H]1O"
    valid, reason = is_oligosaccharide(test_smiles)
    print(valid, reason)