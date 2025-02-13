"""
Classifies: CHEBI:23053 catechin
"""
"""
Classifies: Catechin class – Members of the class of hydroxyflavan that have a flavan-3-ol (catechin) scaffold.
Heuristic:
  • Add explicit hydrogens so that –OH groups are detectable.
  • Search for a six‐membered, non‐aromatic ring that contains exactly one oxygen atom (candidate “C ring”).
  • Check that at least one carbon on that ring is substituted with a hydroxyl (-OH) group.
  • Verify that the candidate ring is fused with an aromatic ring (the “A ring”) sharing at least 2 atoms.
  • Verify that at least one additional aromatic ring (“B ring”) is attached as a substituent (attached to candidate ring atoms not part of the fused system).
Due to diversity in substitution this algorithm is heuristic.
"""

from rdkit import Chem

def is_catechin(smiles: str):
    """
    Determines if a molecule belongs to the catechin (flavan-3-ol) class based on its SMILES string.
    The heuristic searches for a non‐aromatic six‐membered oxygen-containing ring (the C ring)
    fused to an aromatic ring (the A ring) and bearing an additional aromatic substituent (the B ring).
    Explicit hydrogens are added so that hydroxyl substitutions can be identified.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a catechin scaffold, False otherwise.
        str: Reason describing the classification decision.
    """
    # Parse the SMILES string and add explicit hydrogens
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    mol = Chem.AddHs(mol)
    
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()  # each ring is a tuple of atom indices
    
    # Helper: check whether a ring (list of atom indices) is fully aromatic.
    def is_aromatic_ring(ring):
        for idx in ring:
            if not mol.GetAtomWithIdx(idx).GetIsAromatic():
                return False
        return True
    
    # Iterate over all rings looking for a candidate six-membered “C ring”
    for ring in atom_rings:
        if len(ring) != 6:
            continue  # we need a 6-membered ring
        candidate_atoms = [mol.GetAtomWithIdx(i) for i in ring]
        # Count oxygen atoms in the ring (should be exactly one)
        oxy_count = sum(1 for atom in candidate_atoms if atom.GetAtomicNum() == 8)
        if oxy_count != 1:
            continue

        # Ensure the candidate ring is non-aromatic.
        if any(atom.GetIsAromatic() for atom in candidate_atoms):
            continue

        # Check if at least one carbon in the candidate ring has a hydroxyl (-OH) group.
        has_hydroxyl = False
        for atom in candidate_atoms:
            # We focus on carbons (atomic number 6). (The ring oxygen is not used for –OH check.)
            if atom.GetAtomicNum() != 6:
                continue
            # For each neighboring atom, check if it is an -OH group.
            for nbr in atom.GetNeighbors():
                # Skip if neighbor is in the candidate ring
                if nbr.GetIdx() in ring:
                    continue
                # If the substituent is oxygen with at least one hydrogen, count it as hydroxyl.
                if nbr.GetAtomicNum() == 8:
                    # Get explicit hydrogen count
                    # Sometimes the hydrogen count is stored explicitly since we called AddHs.
                    if nbr.GetTotalNumHs() >= 1:
                        has_hydroxyl = True
                        break
            if has_hydroxyl:
                break
        if not has_hydroxyl:
            continue  # candidate ring must have a hydroxyl substituent
        
        # Look for a fused aromatic ring ("A ring"): one that shares at least 2 atoms with the candidate ring.
        fused_aromatic = None
        for other_ring in atom_rings:
            if other_ring == ring:
                continue
            # Check intersection
            common = set(ring) & set(other_ring)
            if len(common) >= 2 and is_aromatic_ring(other_ring):
                fused_aromatic = set(other_ring)
                break
        if fused_aromatic is None:
            continue  # need a fused aromatic ring
        
        # Look for an additional aromatic ring ("B ring") as a substituent.
        # We require that at least one atom of the candidate ring has a neighbor
        # that belongs to an aromatic ring that is NOT the fused aromatic ring.
        b_ring_found = False
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            for nbr in atom.GetNeighbors():
                nbr_idx = nbr.GetIdx()
                # If neighbor is in candidate or in fused aromatic, skip.
                if nbr_idx in ring or nbr_idx in fused_aromatic:
                    continue
                # Check each ring that neighbor belongs to.
                for other_ring in atom_rings:
                    if nbr_idx in other_ring:
                        # Skip rings that are identical to the fused aromatic ring
                        if set(other_ring) == fused_aromatic:
                            continue
                        if is_aromatic_ring(other_ring):
                            b_ring_found = True
                            break
                if b_ring_found:
                    break
            if b_ring_found:
                break
        if not b_ring_found:
            continue  # no additional aromatic substituent found
        
        # If all conditions are met, we assume the molecule contains a catechin (flavan-3-ol) scaffold.
        return True, ("Molecule contains a catechin (flavan-3-ol) scaffold: a non-aromatic six-membered heterocycle "
                      "with one oxygen and a hydroxyl substituent, fused with an aromatic ring and bearing an additional aromatic group.")
    
    return False, "Molecule does not appear to contain the required flavan-3-ol (catechin) scaffold."

# Example usage:
if __name__ == "__main__":
    # Test with one of the provided examples: (-)-catechin
    test_smiles = "O[C@@H]1Cc2c(O)cc(O)cc2O[C@H]1c1ccc(O)c(O)c1"
    decision, reason = is_catechin(test_smiles)
    print("Is catechin?", decision)
    print("Reason:", reason)