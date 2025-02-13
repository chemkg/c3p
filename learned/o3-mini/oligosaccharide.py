"""
Classifies: CHEBI:50699 oligosaccharide
"""
"""
Classifies: Oligosaccharide
Definition: A compound in which monosaccharide units are joined by glycosidic linkages.
Heuristic improvements:
  1. Sugar rings are identified as 5‐membered (furanose) or 6‐membered (pyranose) rings that contain
     exactly one ring oxygen and the remaining ring atoms are mostly carbon. However, one nitrogen is allowed.
  2. At least two such sugar rings must be present and they must be connected by a glycosidic linkage.
  3. The total atoms in these rings should comprise a majority (here ≥50%) of the molecule’s heavy atoms.
Note: This is a heuristic approach and may not cover all edge cases.
"""

from rdkit import Chem

def is_oligosaccharide(smiles: str):
    """
    Determines if a molecule is an oligosaccharide based on its SMILES string.
    
    Approach:
      - Parse the molecule and get its ring information.
      - Identify candidate sugar rings (size 5 or 6) that have exactly one oxygen atom.
        For rings of size 6, we allow either 5 carbons (pyranose) or 4 carbons plus one nitrogen
        (e.g. N-acetyl sugars). Similarly for 5-membered rings.
      - Require that at least two sugar rings are found.
      - Build a connectivity graph among the sugar rings. Two rings are “connected”
        if (a) any atom in one ring is directly bonded to an atom in the other or
        (b) an atom in one ring is bonded to an exocyclic oxygen that in turn is bonded to an atom in the other ring.
      - Check that all candidate sugar rings form one connected subgraph.
      - Finally, compute the fraction of heavy atoms in the molecule that belong to these sugar rings.
        For a valid “oligosaccharide” the sugar portion should dominate (here we require ≥50%).
    
    Args:
       smiles (str): SMILES string of the molecule
       
    Returns:
       bool: True if molecule is an oligosaccharide, False otherwise.
       str: Reason for classification.
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    heavy_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1]
    total_heavy = len(heavy_atoms)
    
    ring_info = mol.GetRingInfo().AtomRings()
    sugar_rings = []  # list of sets of atom indices
    
    # Examine each ring to see if it is a candidate sugar ring.
    for ring in ring_info:
        if len(ring) not in (5, 6):
            continue
        o_count = 0
        c_count = 0
        n_count = 0
        # Count atoms that are part of the ring and what they are
        for idx in ring:
            atomic_num = mol.GetAtomWithIdx(idx).GetAtomicNum()
            if atomic_num == 8:
                o_count += 1
            elif atomic_num == 6:
                c_count += 1
            elif atomic_num == 7:
                n_count += 1
            
        # For a typical sugar ring: should have exactly one oxygen.
        # For a pure pyranose: 6-membered: 1 oxygen and 5 carbons.
        # Allow for one nitrogen to substitute a carbon (e.g. in amino sugars).
        if len(ring) == 6:
            if o_count == 1 and (c_count == 5 or (c_count == 4 and n_count == 1)):
                sugar_rings.append(set(ring))
        elif len(ring) == 5:
            if o_count == 1 and (c_count == 4 or (c_count == 3 and n_count == 1)):
                sugar_rings.append(set(ring))
    
    if len(sugar_rings) < 2:
        return False, f"Fewer than 2 candidate sugar rings detected (found {len(sugar_rings)})."
    
    # Build connectivity graph among sugar rings.
    # Nodes are indices in sugar_rings; an edge exists if two rings are connected by a glycosidic bond.
    n = len(sugar_rings)
    graph = {i: set() for i in range(n)}
    
    for i in range(n):
        for j in range(i+1, n):
            connected = False
            # Check for direct bond between any atom in ring i and any in ring j.
            for bond in mol.GetBonds():
                a_idx = bond.GetBeginAtomIdx()
                b_idx = bond.GetEndAtomIdx()
                if (a_idx in sugar_rings[i] and b_idx in sugar_rings[j]) or (a_idx in sugar_rings[j] and b_idx in sugar_rings[i]):
                    connected = True
                    break
            if not connected:
                # Check for bridging oxygen:
                # For each atom in ring i, check if one of its neighbors (not in ring i) is oxygen,
                # and then if that oxygen is bonded to an atom in ring j.
                for a_idx in sugar_rings[i]:
                    atom_a = mol.GetAtomWithIdx(a_idx)
                    for nbr in atom_a.GetNeighbors():
                        if nbr.GetIdx() in sugar_rings[i]:
                            continue
                        if nbr.GetAtomicNum() != 8:
                            continue
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
    
    # Check if the connectivity graph is connected (all rings in one network).
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
    
    # Finally, check that the sugar portion accounts for a majority of the molecule.
    sugar_atoms = set()
    for ring in sugar_rings:
        sugar_atoms.update(ring)
    fraction = len(sugar_atoms) / total_heavy if total_heavy > 0 else 0
    if fraction < 0.5:
        return False, f"Sugar rings detected but they only constitute {fraction*100:.1f}% of heavy atoms; likely a glycoside rather than an oligosaccharide."
    
    return True, "Molecule has multiple connected sugar rings that constitute the majority of heavy atoms; likely an oligosaccharide."

# Example usage:
if __name__ == '__main__':
    # alpha-mannobiose (a disaccharide) example:
    test_smiles = "OC[C@H]1O[C@@H](O[C@@H]2[C@@H](CO)O[C@H](O)[C@@H](O)[C@H]2O)[C@@H](O)[C@@H](O)[C@@H]1O"
    valid, reason = is_oligosaccharide(test_smiles)
    print(valid, reason)