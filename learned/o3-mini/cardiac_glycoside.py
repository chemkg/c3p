"""
Classifies: CHEBI:83970 cardiac glycoside
"""
"""
Classifies: Cardiac glycoside
Definition: Steroid lactones containing sugar residues that act on the contractile force of the cardiac muscles.

Heuristic approach:
  1. Detect a steroid nucleus by looking for a fused ring system composed of ≥4 rings that share at least 2 atoms.
     When found, store the union of the atom indices of these rings for later connection checks.
  2. Detect a lactone (butenolide) ring by scanning for a 5‐membered ring (ring size 5) that contains at least one oxygen atom and 
     at least one carbon in that ring that is connected via a double bond to an external oxygen (C=O).
  3. Detect a sugar residue by scanning for rings of size 5 or 6 that have exactly one ring oxygen and at least 2 exocyclic –OH groups.
     In addition, at least one atom of the candidate sugar ring must be directly attached (via any bond) to an atom of the steroid nucleus.
     
Note: This is a heuristic approach. It does not guarantee perfect discrimination between cardiac glycosides and other molecules.
"""

from rdkit import Chem

def is_cardiac_glycoside(smiles: str):
    """
    Determines if a molecule qualifies as a cardiac glycoside based on its SMILES string.
    
    Heuristic criteria:
      (1) A steroid nucleus: a set of ≥4 fused rings (fused means sharing at least 2 atoms). 
          The atoms comprising that fused system are stored.
      (2) A lactone ring: a 5-membered ring that contains at least one ring oxygen and at least one carbon 
          (in the ring) that has a double bond to an oxygen outside the ring.
      (3) At least one sugar residue: a 5- or 6-membered ring that (a) has exactly one ring oxygen,
          (b) has at least two exocyclic hydroxyl groups attached to ring atoms and (c) is directly attached 
          to the steroid nucleus.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule meets all criteria, False otherwise.
        str: Explanation of the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    
    # (1) Check for fused steroid nucleus.
    if len(rings) < 4:
        return False, "Fewer than 4 rings found in molecule; no steroid nucleus."
    
    n = len(rings)
    adjacency = {i: set() for i in range(n)}
    # Two rings are considered fused if they share at least 2 atoms.
    for i in range(n):
        for j in range(i+1, n):
            if len(set(rings[i]).intersection(rings[j])) >= 2:
                adjacency[i].add(j)
                adjacency[j].add(i)
    
    visited = set()
    steroid_nucleus_found = False
    steroid_atoms = set()
    def dfs(node, comp):
        comp.add(node)
        for neighbor in adjacency[node]:
            if neighbor not in comp:
                dfs(neighbor, comp)
        return comp
    
    for i in range(n):
        if i not in visited:
            comp = dfs(i, set())
            # if fused component has at least 4 rings then mark as steroid nucleus
            if len(comp) >= 4:
                steroid_nucleus_found = True
                # Collect all atoms in these rings for later sugar attachment check.
                for ring_idx in comp:
                    steroid_atoms.update(rings[ring_idx])
                break
            visited.update(comp)
    if not steroid_nucleus_found:
        return False, "No fused steroid nucleus (≥4 fused rings) detected"
    
    # (2) Lactone (butenolide) detection: look for a 5-membered ring that has at least one ring oxygen and a carbon
    # in the ring that has an external C=O bond.
    lactone_found = False
    for ring in rings:
        if len(ring) == 5:
            ring_oxygens = [idx for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8]
            if len(ring_oxygens) < 1:
                continue
            # Look for a carbon in ring that is bonded (double bond) to an oxygen not in the ring.
            for idx in ring:
                atom = mol.GetAtomWithIdx(idx)
                if atom.GetAtomicNum() == 6:  # carbon
                    for bond in atom.GetBonds():
                        # Check bond double order (approximate by GetBondTypeAsDouble)
                        if bond.GetBondTypeAsDouble() == 2.0:
                            nbr = bond.GetOtherAtom(atom)
                            if nbr.GetAtomicNum() == 8 and nbr.GetIdx() not in ring:
                                lactone_found = True
                                break
                    if lactone_found:
                        break
            if lactone_found:
                break
    if not lactone_found:
        return False, "No lactone (butenolide) ring detected"
    
    # (3) Sugar residue detection.
    sugar_found = False
    for ring in rings:
        if len(ring) not in (5, 6):
            continue
        ring_atoms = [mol.GetAtomWithIdx(i) for i in ring]
        # For typical sugars (pyranoses/furanoses) we expect exactly one ring oxygen.
        ring_oxygen_count = sum(1 for atom in ring_atoms if atom.GetAtomicNum() == 8)
        if ring_oxygen_count != 1:
            continue
        
        # Count exocyclic hydroxyl substituents on the ring atoms.
        hydroxyl_count = 0
        for atom in ring_atoms:
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() not in ring and nbr.GetAtomicNum() == 8:
                    # Consider it a hydroxyl if the oxygen has at least one hydrogen.
                    if nbr.GetTotalNumHs() > 0:
                        hydroxyl_count += 1
        # Relax threshold to at least 2 hydroxyl groups.
        if hydroxyl_count < 2:
            continue
        
        # Additionally, check that the candidate sugar ring is attached to the steroid nucleus.
        attached_to_steroid = False
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() not in ring and nbr.GetIdx() in steroid_atoms:
                    attached_to_steroid = True
                    break
            if attached_to_steroid:
                break
        if attached_to_steroid:
            sugar_found = True
            break
    if not sugar_found:
        return False, "No sugar residue (glycosidic moiety) detected that is directly attached to the steroid nucleus"
    
    msg = ("Detected a steroid nucleus (≥4 fused rings), a lactone (butenolide) ring, and at least one sugar residue that is attached "
           "to the nucleus – criteria for a cardiac glycoside structure.")
    return True, msg


# Example usage:
if __name__ == "__main__":
    # Test with one example of a known cardiac glycoside (Erychroside)
    test_smiles = "O[C@@]12[C@]3([C@@]([C@@]4([C@](O)(CC3)C[C@@H](O[C@@H]5O[C@@H]([C@@H](O[C@@H]6OC[C@@H](O)[C@H](O)[C@H]6O)[C@@H](O)C5)C)CC4)C=O)(CC[C@@]1([C@H](CC2)C=7COC(=O)C7)C)[H])[H]"
    result, reason = is_cardiac_glycoside(test_smiles)
    print("Is cardiac glycoside?", result)
    print("Reason:", reason)