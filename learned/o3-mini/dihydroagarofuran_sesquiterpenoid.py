"""
Classifies: CHEBI:71548 dihydroagarofuran sesquiterpenoid
"""
#!/usr/bin/env python
"""
Classifies: dihydroagarofuran sesquiterpenoid
Definition: Any sesquiterpenoid with a dihydroagarofuran skeleton.

Improved heuristic:
  1. Identify candidate tetrahydrofuran rings (5-membered rings with exactly 1 oxygen and 4 carbons).
  2. For each candidate, grow a fused ring system by including any ring that shares at least one atom.
  3. Count the number of carbon atoms in the union of the atoms in the fused system.
  4. If the carbon count is in the range consistent with a sesquiterpenoid core (here ~14-19 C atoms),
     then classify the molecule as a dihydroagarofuran sesquiterpenoid.
"""

from rdkit import Chem

def is_dihydroagarofuran_sesquiterpenoid(smiles: str):
    """
    Determines if a molecule is a dihydroagarofuran sesquiterpenoid based on its SMILES string.

    This improved function uses a relaxed version of the previous heuristic:
      - It searches for a tetrahydrofuran ring (defined as a 5-membered ring with 1 oxygen and 4 carbons).
      - For each candidate THF ring, it collects all rings that are fused (by sharing at least one atom)
        into a connected set.
      - It counts the carbons in the union of all atoms present in this fused ring system.
      - If the total carbon count falls in the range (~14 to 19),
        the molecule is classified as a dihydroagarofuran sesquiterpenoid.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is classified as a dihydroagarofuran sesquiterpenoid, False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse the molecule from SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get ring information from the molecule
    ring_info = mol.GetRingInfo()
    if not ring_info or not ring_info.AtomRings():
        return False, "No rings detected in the molecule"
    rings = ring_info.AtomRings()  # list of tuples of atom indices

    # Identify candidate tetrahydrofuran (THF) rings: 5-membered ring with exactly one oxygen and four carbons.
    candidate_thf = []
    for ring in rings:
        if len(ring) == 5:
            oxygen_count = 0
            carbon_count = 0
            for idx in ring:
                atom_sym = mol.GetAtomWithIdx(idx).GetSymbol()
                if atom_sym == 'O':
                    oxygen_count += 1
                elif atom_sym == 'C':
                    carbon_count += 1
            if oxygen_count == 1 and carbon_count == 4:
                candidate_thf.append(set(ring))
    
    if not candidate_thf:
        return False, "No tetrahydrofuran (5-membered, 1O/4C) ring found"
    
    # For each candidate THF ring, build the fused ring system using a relaxed connectivity rule:
    # Any ring sharing at least one atom with the candidate or with a ring already in the
    # fused set will be included.
    for candidate in candidate_thf:
        fused_rings = [candidate]
        updated = True
        # Loop until no new ring is added
        while updated:
            updated = False
            for ring in rings:
                ring_set = set(ring)
                # Skip rings already in the fused system (using set equality)
                if any(ring_set == existing for existing in fused_rings):
                    continue
                # If the ring shares at least one atom with any ring already in the fused system, add it.
                if any(len(ring_set.intersection(existing)) >= 1 for existing in fused_rings):
                    fused_rings.append(ring_set)
                    updated = True
        
        # Compute the union of all atoms in the fused ring system.
        fused_atom_union = set()
        for ring_set in fused_rings:
            fused_atom_union = fused_atom_union.union(ring_set)
        
        # Count the number of carbons in this fused ring system.
        core_carbons = sum(1 for idx in fused_atom_union if mol.GetAtomWithIdx(idx).GetSymbol() == 'C')
        
        # We expect a dihydroagarofuran sesquiterpenoid core to have about 15 carbons.
        # By testing, a range of 14 to 19 carbons seems to capture the true positives
        # while rejecting many false positives.
        if 14 <= core_carbons <= 19:
            reason = ("Fused ring system detected with a tetrahydrofuran ring and ~{} core carbons "
                      "(expected for a sesquiterpenoid).".format(core_carbons))
            return True, reason
            
    return False, "No fused ring system consistent with a dihydroagarofuran skeleton was found"

# If run as a script, execute one example test.
if __name__ == "__main__":
    # Example: one provided SMILES for wilforidine (a known dihydroagarofuran sesquiterpenoid).
    test_smiles = "[C@@]12([C@@H]([C@@H]([C@@]3([C@H]([C@]14[C@]([C@H]([C@@H]([C@@H]2OC(=O)C)O)OC(C(CCC5=NC=CC=C5C(OC[C@@]3(O4)C)=O)(C)O)=O)(O)C)OC(C)=O)[H])OC(C)=O)OC(=O)C)COC(C)=O"
    result, explanation = is_dihydroagarofuran_sesquiterpenoid(test_smiles)
    print("Result:", result)
    print("Explanation:", explanation)