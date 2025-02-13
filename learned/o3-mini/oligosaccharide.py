"""
Classifies: CHEBI:50699 oligosaccharide
"""
"""
Classifies: Oligosaccharide
Definition: A compound in which monosaccharide units are joined by glycosidic linkages.
Heuristic improvements in this version:
  1. Candidate sugar rings are identified as either 5‐membered (furanose) or 6‐membered (pyranose) rings that are fully saturated 
     (all bonds in the ring are single) and contain exactly one ring oxygen (with allowance for one nitrogen instead of a carbon).
  2. For each candidate ring, at least two ring carbons must have an exocyclic hydroxyl group 
     (an oxygen neighbor having at least one attached hydrogen).
  3. At least two candidate sugar rings must be found and they must be connected via what looks like a glycosidic linkage.
     Two rings are considered “connected” only if (a) a direct bond joins atoms from the two rings and one of the atoms is oxygen,
     or (b) an oxygen atom (with degree exactly 2) bridges a carbon in one ring to a carbon in the other.
  4. Finally, the atoms in these rings should account for at least 60% of the molecule’s heavy atoms.
Note: This is still a heuristic approach that will not cover all edge cases.
"""

from rdkit import Chem

def is_oligosaccharide(smiles: str):
    """
    Determines if a molecule is an oligosaccharide based on its SMILES string.
    Uses several heuristics:
      - Identify candidate sugar rings (5- or 6-membered) that are fully saturated and have exactly one ring oxygen.
      - Each candidate ring must have at least 2 carbon atoms with exocyclic -OH groups.
      - At least 2 candidate sugar rings must be present and connected by glycosidic linkages.
        Here a glycosidic connection is accepted only if:
          (a) a bond joins an atom from one ring to an atom in the other and at least one partner is oxygen (and the bond is single), or 
          (b) a bridging oxygen (with degree exactly 2) bonds to a carbon in one ring and a carbon in the other.
      - The sugar rings overall must account for at least 60% of the heavy (non-hydrogen) atoms.
      
    Args:
       smiles (str): SMILES string of the molecule.
       
    Returns:
       bool: True if the molecule is classified as an oligosaccharide, False otherwise.
       str: A reason for the classification decision.
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Count heavy atoms (atomic number > 1)
    heavy_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1]
    total_heavy = len(heavy_atoms)
    
    ring_info = mol.GetRingInfo().AtomRings()
    sugar_rings = []  # each candidate stored as a set of atom indices
    
    # Loop over each ring in the molecule
    for ring in ring_info:
        if len(ring) not in (5, 6):
            continue
        ring_set = set(ring)
        
        # Check that every bond within the ring is a single bond.
        saturated = True
        for bond in mol.GetBonds():
            a_idx = bond.GetBeginAtomIdx()
            b_idx = bond.GetEndAtomIdx()
            if a_idx in ring_set and b_idx in ring_set:
                if bond.GetBondType() != Chem.BondType.SINGLE:
                    saturated = False
                    break
        if not saturated:
            continue
        
        # Count atoms by atomic number in the ring
        o_count = c_count = n_count = 0
        for idx in ring:
            atomic_num = mol.GetAtomWithIdx(idx).GetAtomicNum()
            if atomic_num == 8:
                o_count += 1
            elif atomic_num == 6:
                c_count += 1
            elif atomic_num == 7:
                n_count += 1
        
        # Typical sugar ring: exactly one oxygen.
        if len(ring) == 6:
            if not (o_count == 1 and (c_count == 5 or (c_count == 4 and n_count == 1))):
                continue
        elif len(ring) == 5:
            if not (o_count == 1 and (c_count == 4 or (c_count == 3 and n_count == 1))):
                continue
        
        # Check that at least 2 ring carbons show exocyclic hydroxyls.
        hydroxyl_carbons = 0
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() != 6:  # only consider carbons
                continue
            # Look at neighbors that are NOT in the ring
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring_set:
                    continue
                if nbr.GetAtomicNum() == 8 and nbr.GetTotalNumHs() > 0:
                    hydroxyl_carbons += 1
                    break  # count once per carbon
        if hydroxyl_carbons < 2:
            continue

        # Passed all checks: add candidate ring.
        sugar_rings.append(ring_set)
    
    if len(sugar_rings) < 2:
        return False, f"Fewer than 2 candidate sugar rings detected (found {len(sugar_rings)})."
    
    # Define a helper function to check glycosidic connectivity between two candidate rings.
    def rings_connected(ring1, ring2, mol_obj):
        # Option (a): Direct bond between atoms of different rings where at least one atom is oxygen.
        for bond in mol_obj.GetBonds():
            a_idx = bond.GetBeginAtomIdx()
            b_idx = bond.GetEndAtomIdx()
            if (a_idx in ring1 and b_idx in ring2) or (a_idx in ring2 and b_idx in ring1):
                a = mol_obj.GetAtomWithIdx(a_idx)
                b = mol_obj.GetAtomWithIdx(b_idx)
                if (a.GetAtomicNum() == 8 or b.GetAtomicNum() == 8) and bond.GetBondType() == Chem.BondType.SINGLE:
                    return True
        # Option (b): Bridging oxygen.
        for a_idx in ring1:
            a = mol_obj.GetAtomWithIdx(a_idx)
            for nbr in a.GetNeighbors():
                # Look only for oxygen atoms not in ring1
                if nbr.GetAtomicNum() != 8 or nbr.GetIdx() in ring1:
                    continue
                if nbr.GetDegree() != 2:
                    continue
                # For the bridging oxygen, check its neighbors.
                nbr_idxs = [n.GetIdx() for n in nbr.GetNeighbors()]
                # Exclude the current atom a_idx.
                for n_idx in nbr_idxs:
                    if n_idx in ring2:
                        return True
        return False
    
    # Build connectivity graph among the candidate sugar rings.
    n_rings = len(sugar_rings)
    graph = {i: set() for i in range(n_rings)}
    for i in range(n_rings):
        for j in range(i+1, n_rings):
            if rings_connected(sugar_rings[i], sugar_rings[j], mol):
                graph[i].add(j)
                graph[j].add(i)
    
    # Check if all candidate rings are connected (i.e. form one connected component)
    visited = set()
    stack = [0]
    while stack:
        node = stack.pop()
        if node in visited:
            continue
        visited.add(node)
        for nb in graph[node]:
            if nb not in visited:
                stack.append(nb)
    if len(visited) != n_rings:
        return False, "Candidate sugar rings are detected but they are not all connected via glycosidic linkages."
    
    # Count total atoms that are found in any candidate ring.
    ring_atom_indices = set()
    for ring in sugar_rings:
        ring_atom_indices.update(ring)
    fraction = len(ring_atom_indices) / total_heavy if total_heavy else 0
    if fraction < 0.6:
        return False, f"Sugar rings detected but they constitute only {fraction*100:.1f}% of heavy atoms; likely a glycoside rather than an oligosaccharide."
    
    return True, ("Molecule has multiple connected, saturated sugar rings with exocyclic hydroxyl groups " 
                  "and these rings account for the majority (≥60%) of heavy atoms; likely an oligosaccharide.")

# Example usage:
if __name__ == '__main__':
    # Example: alpha-mannobiose (a disaccharide)
    test_smiles = "OC[C@H]1O[C@@H](O[C@@H]2[C@@H](CO)O[C@H](O)[C@@H](O)[C@H]2O)[C@@H](O)[C@@H](O)[C@@H]1O"
    valid, reason = is_oligosaccharide(test_smiles)
    print(valid, reason)