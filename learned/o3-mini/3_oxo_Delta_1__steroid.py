"""
Classifies: CHEBI:20156 3-oxo-Delta(1) steroid
"""
"""
Classifies: CHEBI:3-oxo-Δ(1) steroid
Definition: Any 3-oxo steroid that contains a double bond between positions 1 and 2.
The improved heuristic:
  - Require that the molecule has at least 4 rings.
  - Identify six-membered rings that are fused with at least one other 5- or 6-membered ring.
  - In each such candidate ring, look for a ketone carbon (C=O) that is not an aldehyde.
  - Check that the same ring contains at least one non-carbonyl C=C double bond, serving as a proxy 
    for the Δ(1) double bond.
Note: This heuristic does not perfectly capture all edge cases.
"""

from rdkit import Chem

def is_3_oxo_Delta_1__steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-Δ(1) steroid based on its SMILES string.

    Heuristic criteria:
      1. The molecule must have at least 4 rings (the steroid nucleus is typically tetracyclic).
      2. Look for six-membered rings that are fused to at least one other 5- or 6-membered ring 
         (i.e. share 2 or more atoms).
      3. In such a six-membered ring, search for a ketone carbon (a carbon double bonded to oxygen)
         that is coupled to at least two other carbon atoms (to avoid aldehydes).
      4. Also in that ring, check for at least one non-carbonyl C=C double bond (proxy for the Δ(1) bond).
    
    Args:
       smiles (str): SMILES string representing the molecule.
    
    Returns:
       bool: True if it qualifies as a 3-oxo-Δ(1) steroid, else False.
       str: Explanation for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get ring information from the molecule.
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    if len(rings) < 4:
        return False, "Molecule does not have at least 4 rings required for a steroid nucleus"

    # We'll focus on candidate rings that are exactly six-membered.
    candidate_rings = [r for r in rings if len(r) == 6]
    if not candidate_rings:
        return False, "No six-membered rings found"

    # Define a helper function to test if a ring is fused with another ring.
    def is_fused(ring):
        ring_set = set(ring)
        # Look for a different ring (of size 5 or 6) that shares at least 2 atoms.
        for other in rings:
            if other == ring:
                continue
            if len(other) not in (5, 6):
                continue
            if len(ring_set.intersection(other)) >= 2:
                return True
        return False

    # For each candidate six-membered ring that is fused,
    # look for a ketone carbon and a C=C double bond (non-carbonyl) in that ring.
    for ring in candidate_rings:
        if not is_fused(ring):
            continue  # skip rings that are isolated
        ring_set = set(ring)

        candidate_ketone_found = False
        # iterate over atoms in the ring to find a ketone (C=O) center.
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetSymbol() != "C":
                continue
            # Look for a double bond from this carbon to an oxygen.
            for bond in atom.GetBonds():
                # Check for a double bond.
                if bond.GetBondTypeAsDouble() == 2:
                    other = bond.GetOtherAtom(atom)
                    if other.GetSymbol() == "O":
                        # Ensure that the carbon is bonded to at least two carbons (ketone vs aldehyde check).
                        carbon_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 6]
                        if len(carbon_neighbors) >= 2:
                            candidate_ketone_found = True
                            break
            if candidate_ketone_found:
                break

        # If no ketone candidate in this ring, go on to the next ring.
        if not candidate_ketone_found:
            continue

        # Now, check this candidate ring for a non-carbonyl C=C double bond.
        found_delta = False
        # iterate over bonds in the molecule and focus on bonds where both atoms lie in this ring.
        for bond in mol.GetBonds():
            a1_idx = bond.GetBeginAtomIdx()
            a2_idx = bond.GetEndAtomIdx()
            if a1_idx in ring_set and a2_idx in ring_set:
                # Must be a double bond.
                if bond.GetBondTypeAsDouble() == 2:
                    # Ensure it is a bond between two carbons.
                    a1 = mol.GetAtomWithIdx(a1_idx)
                    a2 = mol.GetAtomWithIdx(a2_idx)
                    if a1.GetSymbol() == "C" and a2.GetSymbol() == "C":
                        # Do not consider the carbonyl bond (C=O) because that was used for the ketone check.
                        # (Since the oxygen isn't in the ring, this bond will not be counted.)
                        found_delta = True
                        break
        if found_delta:
            return True, ("Molecule contains a fused steroid nucleus (>=4 rings) with a six-membered ring "
                          "bearing a ketone (3-oxo) and a non-carbonyl C=C double bond (proxy for Δ(1) bond).")
    return False, "No candidate 3-oxo-Δ(1) steroid nucleus found"


# Test cases (you can remove or comment these out)
if __name__ == "__main__":
    test_smiles = [
        # paraminabeolide B (expected True)
        "C[C@@H]([C@H]1CC[C@H]2[C@@H]3CCC4=CC(=O)C=C[C@]4(C)[C@H]3CC[C@]12COC(C)=O",
        # helvolic acid methyl ester (expected True)
        "C=1[C@@]2([C@@]3(CC[C@@]/4([C@@]([C@]3(C([C@H]([C@]2([C@@H](C(C1)=O)C)[H])OC(=O)C)=O)C)(C[C@@H](\\C4=C(\\CCC=C(C)C)/C(=O)OC)OC(=O)C)C)[H])[H])C",
        # estra-1,5(10)-diene-3,4,17-trione (expected True)
        "C[C@]12CC[C@H]3[C@@H](CCC4=C3C=CC(=O)C4=O)[C@@H]1CCC2=O",
    ]
    for smi in test_smiles:
        res, reason = is_3_oxo_Delta_1__steroid(smi)
        print("SMILES:", smi)
        print("Classification:", res)
        print("Explanation:", reason)
        print("----------")