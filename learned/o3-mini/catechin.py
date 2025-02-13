"""
Classifies: CHEBI:23053 catechin
"""
"""
Classifies: Catechin class – Members of the class of hydroxyflavan that have a flavan-3-ol (catechin) scaffold.
Heuristic:
  • Look for a non‐aromatic six‐membered ring in the molecule that contains exactly one oxygen atom.
  • Check if one or more ring carbons in that ring have a hydroxyl (–OH) substituent.
  • Verify that the ring is fused with an aromatic ring (sharing at least 2 atoms – the “A ring”).
  • Verify that at least one additional aromatic ring (a “B ring”) is attached as a substituent.
Due to diversity in substitution, this algorithm is a heuristic.
"""

from rdkit import Chem

def is_catechin(smiles: str):
    """
    Determines if a molecule belongs to the catechin (flavan-3-ol) class based on its SMILES string.
    It does so by searching for a six‐membered, non‐aromatic (saturated) oxygen‐containing ring 
    (the flavan “C ring”) that:
      - contains exactly one oxygen,
      - carries at least one hydroxyl (–OH) group on one of its carbons,
      - is fused (shares at least 2 atoms) with an aromatic ring (the “A ring”),
      - has at least one additional attached aromatic ring (the “B ring”).

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is classified as a catechin (flavan-3-ol), False otherwise.
        str: Reason describing the classification decision.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."

    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()  # tuple of tuples, each is a ring (list of atom indices)

    # Helper: test if a ring (given as indices) is entirely aromatic.
    def is_aromatic_ring(ring):
        # A ring is aromatic if every atom in it is flagged as aromatic.
        for idx in ring:
            if not mol.GetAtomWithIdx(idx).GetIsAromatic():
                return False
        return True

    # Search among all rings for a candidate "C-ring": six-membered, non-aromatic, exactly one oxygen.
    for ring in atom_rings:
        if len(ring) != 6:
            continue  # not six-membered
        candidate_atoms = [mol.GetAtomWithIdx(i) for i in ring]
        # Count the number of oxygen atoms in the candidate ring
        oxy_count = sum(1 for atom in candidate_atoms if atom.GetAtomicNum() == 8)
        if oxy_count != 1:
            continue  # must contain exactly one oxygen (the heteroatom)
        # In a catechin's C-ring, the carbons should be non-aromatic.
        if any(atom.GetAtomicNum() == 6 and atom.GetIsAromatic() for atom in candidate_atoms):
            continue

        # Check that at least one carbon in the candidate ring bears a hydroxyl group.
        hydroxyl_found = False
        for atom in candidate_atoms:
            if atom.GetAtomicNum() == 6:
                # Look at neighbors not in the ring. A hydroxyl group is an -OH: oxygen with at least one hydrogen.
                for nbr in atom.GetNeighbors():
                    if nbr.GetAtomicNum() == 8 and nbr.GetTotalNumHs() > 0:
                        hydroxyl_found = True
                        break
            if hydroxyl_found:
                break
        if not hydroxyl_found:
            continue

        # Look for a fused aromatic ring (“A ring”) that shares at least 2 atoms with the candidate ring.
        fused_aromatic = None
        for other_ring in atom_rings:
            if other_ring == ring:
                continue
            # Check if the intersection contains at least 2 atoms.
            if len(set(ring) & set(other_ring)) >= 2:
                # Consider this other ring aromatic if all its atoms are aromatic.
                if is_aromatic_ring(other_ring):
                    fused_aromatic = other_ring
                    break
        if fused_aromatic is None:
            continue  # no fused aromatic ring found

        # Look for an additional aromatic substituent (“B ring”) attached to the candidate ring.
        # We require that at least one atom of the candidate ring has a neighbor that is part of an aromatic ring
        # that is not the fused aromatic ring.
        b_ring_found = False
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            for nbr in atom.GetNeighbors():
                nbr_idx = nbr.GetIdx()
                if nbr_idx in ring:
                    continue  # neighbor is within the candidate ring
                # Check all rings that nbr participates in.
                for other_ring in atom_rings:
                    if nbr_idx in other_ring:
                        # Skip the fused aromatic ring we already found
                        if fused_aromatic is not None and set(other_ring) == set(fused_aromatic):
                            continue
                        # Accept if other_ring is aromatic
                        if is_aromatic_ring(other_ring):
                            b_ring_found = True
                            break
                if b_ring_found:
                    break
            if b_ring_found:
                break
        if not b_ring_found:
            continue  # did not find an extra aromatic substituent

        # If all conditions are met, we assume the molecule contains a catechin (flavan-3-ol) scaffold.
        return True, ("Molecule contains a catechin (flavan-3-ol) scaffold: a non-aromatic six-membered heterocycle "
                      "with one oxygen and a hydroxyl group, fused to an aromatic ring and bearing an extra aromatic substituent.")

    return False, "Molecule does not appear to contain the required flavan-3-ol (catechin) scaffold."

# Example usage:
if __name__ == "__main__":
    # Test with one of the provided examples: (–)-catechin
    test_smiles = "O[C@@H]1Cc2c(O)cc(O)cc2O[C@H]1c1ccc(O)c(O)c1"
    decision, reason = is_catechin(test_smiles)
    print("Is catechin?", decision)
    print("Reason:", reason)