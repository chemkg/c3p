"""
Classifies: CHEBI:83970 cardiac glycoside
"""
"""
Classifies: Cardiac glycoside
Definition: Steroid lactones containing sugar residues that act on the contractile force of the cardiac muscles.

Heuristic approach:
  1. Detect a steroid nucleus by checking for a fused ring system of 4 or more rings.
  2. Detect a lactone (butenolide) ring by scanning 5-membered rings for a ring oxygen and 
     a carbon (in the ring) that bears an external carbonyl (C=O) double bond.
  3. Detect a sugar residue by scanning rings of size 5 or 6 that have one ring oxygen and at least
     two hydroxyl (-OH) substituents on the ring atoms.
     
Note: This is a heuristic approach. The code does not guarantee perfect accuracy.
"""

from rdkit import Chem

def is_cardiac_glycoside(smiles: str):
    """
    Determines if a molecule is a cardiac glycoside based on its SMILES string.
    Cardiac glycosides are defined as steroid lactones containing sugar residues.
    
    Heuristic criteria:
      - Steroid nucleus: fused ring system with at least 4 rings that share bonds.
      - Lactone ring: a 5-membered ring that contains one oxygen and a carbon bearing an external carbonyl.
      - Sugar residue: at least one ring of size 5 or 6 containing one ring oxygen and having at least two hydroxyl substituents.
      
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule qualifies as a cardiac glycoside, False otherwise.
        str: Reason for the classification.
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # ---------------------------
    # (1) Check for fused steroid nucleus.
    # We define a steroid nucleus as a set of at least 4 rings that are fused (share at least 2 atoms).
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
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
    def dfs(node, comp):
        comp.add(node)
        for neighbor in adjacency[node]:
            if neighbor not in comp:
                dfs(neighbor, comp)
        return comp

    steroid_nucleus_found = False
    for i in range(n):
        if i not in visited:
            comp = dfs(i, set())
            if len(comp) >= 4:
                steroid_nucleus_found = True
                break
            visited.update(comp)
    if not steroid_nucleus_found:
        return False, "No fused steroid nucleus (≥4 fused rings) detected"

    # ---------------------------
    # (2) Check for lactone (butenolide) ring.
    # We search for 5-membered rings that (a) contain at least one oxygen in the ring,
    # and (b) have at least one ring carbon that bears a carbonyl double bond (C=O) outside the ring.
    lactone_found = False
    for ring in rings:
        if len(ring) == 5:  # look at 5-membered rings
            # Count oxygen atoms in the ring:
            ring_oxygens = [idx for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8]
            if len(ring_oxygens) < 1:
                continue  # must at least have one oxygen in the ring
            # Now, check every atom in the ring that is a carbon to see if it is attached (via a double bond)
            # to an oxygen that is not part of the ring.
            for idx in ring:
                atom = mol.GetAtomWithIdx(idx)
                if atom.GetAtomicNum() == 6:  # carbon atom
                    for bond in atom.GetBonds():
                        # Look for a double bond
                        if bond.GetBondTypeAsDouble() == 2.0:
                            nbr = bond.GetOtherAtom(atom)
                            # The neighbor should be oxygen and should lie outside the ring.
                            if nbr.GetAtomicNum() == 8 and nbr.GetIdx() not in ring:
                                lactone_found = True
                                break
                    if lactone_found:
                        break
            if lactone_found:
                break
    if not lactone_found:
        return False, "No lactone (butenolide) ring detected"
    
    # ---------------------------
    # (3) Check for sugar residues.
    # Instead of fixed SMARTS, we now look over rings that are 5 or 6 members in size.
    # A typical sugar (pyranose or furanose) ring has one oxygen atom in the ring and several hydroxyl groups attached.
    sugar_found = False
    for ring in rings:
        if len(ring) not in (5, 6):
            continue
        # Count oxygen atoms within the ring:
        ring_atoms = [mol.GetAtomWithIdx(i) for i in ring]
        ring_oxygen_count = sum(1 for atom in ring_atoms if atom.GetAtomicNum() == 8)
        # For a sugar ring, we expect exactly one ring oxygen.
        if ring_oxygen_count != 1:
            continue
        # Now, count the number of hydroxyl substituents on the ring atoms.
        # We consider a substituent hydroxyl if an atom not in the ring is oxygen and has at least one hydrogen.
        hydroxyl_count = 0
        for atom in ring_atoms:
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() not in ring and nbr.GetAtomicNum() == 8:
                    # Use GetTotalNumHs which returns the (implicit + explicit) hydrogen count.
                    if nbr.GetTotalNumHs() > 0:
                        hydroxyl_count += 1
        # We set a loose threshold that a sugar ring will have at least 2 hydroxyl groups.
        if hydroxyl_count >= 2:
            sugar_found = True
            break
    if not sugar_found:
        return False, "No sugar residue (glycosidic moiety) detected"
    
    # If all criteria are met then it is classified as a cardiac glycoside.
    msg = ("Detected a steroid nucleus (≥4 fused rings), a lactone (butenolide) ring, and "
           "at least one sugar residue consistent with a cardiac glycoside structure.")
    return True, msg


# Example usage:
if __name__ == "__main__":
    # Test with one of the provided SMILES (Erychroside).
    test_smiles = "O[C@@]12[C@]3([C@@]([C@@]4([C@](O)(CC3)C[C@@H](O[C@@H]5O[C@@H]([C@@H](O[C@@H]6OC[C@@H](O)[C@H](O)[C@H]6O)[C@@H](O)C5)C)CC4)C=O)(CC[C@@]1([C@H](CC2)C=7COC(=O)C7)C)[H])[H]"
    result, reason = is_cardiac_glycoside(test_smiles)
    print("Is cardiac glycoside?", result)
    print("Reason:", reason)