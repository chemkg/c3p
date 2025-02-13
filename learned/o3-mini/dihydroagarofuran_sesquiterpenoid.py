"""
Classifies: CHEBI:71548 dihydroagarofuran sesquiterpenoid
"""
#!/usr/bin/env python
"""
Classifies: dihydroagarofuran sesquiterpenoid
Definition: Any sesquiterpenoid with a dihydroagarofuran skeleton.
Heuristic (improved):
  1. Identify candidate tetrahydrofuran (THF) rings: 5-membered rings with exactly 1 oxygen and 4 carbon atoms.
  2. For each candidate THF ring, build a fused ring system by grouping rings that share at least two atoms (i.e. a bond).
  3. In that fused system:
       - Exactly one ring should meet the strict THF criteria.
       - The remaining rings must be fully carbocyclic.
       - The total number of rings is 3 or 4.
       - The union of atoms in all fused rings (the “core”) should have 14–16 carbon atoms.
       - The core should contain exactly one oxygen.
  4. If these conditions are met, we classify the molecule as a dihydroagarofuran sesquiterpenoid.
Note: This heuristic is approximate.
"""

from rdkit import Chem

def is_dihydroagarofuran_sesquiterpenoid(smiles: str):
    """
    Determines if a molecule is a dihydroagarofuran sesquiterpenoid based on its SMILES string.

    The function uses an improved heuristic:
      - It looks for candidate tetrahydrofuran (THF) rings (5-membered rings with 1 oxygen and 4 carbons).
      - For each candidate, it gathers rings fused via a common bond (sharing at least 2 atoms).
      - It then requires that in that fused system there is:
            * exactly one THF ring,
            * 2–3 additional fully carbocyclic rings,
            * a total of 3 or 4 rings,
            * and that the "core" (union of fused ring atoms) contains 14–16 carbon atoms and exactly one oxygen.
      - If these conditions are met, the molecule is classified as a dihydroagarofuran sesquiterpenoid.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is classified as a dihydroagarofuran sesquiterpenoid, False otherwise.
        str: Explanation for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    if not atom_rings:
        return False, "No rings detected in the molecule"

    # Helper: Check if a given ring (list of atom indices) fits the strict THF criteria.
    def is_THF(ring):
        if len(ring) != 5:
            return False
        oxygen_count = 0
        carbon_count = 0
        for idx in ring:
            sym = mol.GetAtomWithIdx(idx).GetSymbol()
            if sym == "O":
                oxygen_count += 1
            elif sym == "C":
                carbon_count += 1
            else:
                # If any atom is not C or O, do not consider it a THF candidate.
                return False
        return (oxygen_count == 1 and carbon_count == 4)
    
    # Identify candidate THF rings (using set for easier later handling)
    candidate_thf_sets = []
    for ring in atom_rings:
        if is_THF(list(ring)):
            candidate_thf_sets.append(set(ring))
    
    if not candidate_thf_sets:
        return False, "No tetrahydrofuran (5-membered, 1O/4C) ring found"
    
    # To decide if two rings are fused, we require them to share at least 2 atoms (a bond)
    def rings_are_fused(ring_set_a, ring_set_b):
        return len(ring_set_a.intersection(ring_set_b)) >= 2

    # For a given starting candidate (a candidate THF set), build the fused ring system.
    # We only consider rings fused if they share two or more atoms.
    def build_fused_system(start_ring_set):
        fused_sets = [start_ring_set]
        added = True
        while added:
            added = False
            # iterate over all rings in the molecule
            for ring in atom_rings:
                ring_set = set(ring)
                # Skip rings already in our fused system (exact set match)
                if any(ring_set == f for f in fused_sets):
                    continue
                # If ring_set is fused (shares at least 2 atoms) with any ring already in the system, add it.
                if any(rings_are_fused(ring_set, f) for f in fused_sets):
                    fused_sets.append(ring_set)
                    added = True
        return fused_sets

    # Now analyze each candidate THF's fused system.
    for candidate in candidate_thf_sets:
        fused_rings = build_fused_system(candidate)
        # Count ring types within the fused system.
        thf_count = 0
        carbocyclic_count = 0
        for ring_set in fused_rings:
            ring_list = list(ring_set)
            if is_THF(ring_list):
                thf_count += 1
            else:
                # Fully carbocyclic ring: every atom in the ring must be carbon.
                if all(mol.GetAtomWithIdx(idx).GetSymbol() == "C" for idx in ring_set):
                    carbocyclic_count += 1
        total_rings = len(fused_rings)

        # Get the union of all atom indices from the fused rings.
        fused_atom_union = set()
        for ring_set in fused_rings:
            fused_atom_union |= ring_set

        # Count carbons and oxygens in this union (the core skeleton)
        core_carbons = sum(1 for idx in fused_atom_union if mol.GetAtomWithIdx(idx).GetSymbol() == "C")
        core_oxygens = sum(1 for idx in fused_atom_union if mol.GetAtomWithIdx(idx).GetSymbol() == "O")

        # Our revised conditions:
        #   • exactly one THF ring,
        #   • 2 to 3 additional rings (fully carbocyclic),
        #   • a total fused ring system of 3 or 4 rings,
        #   • core carbons between 14 and 16,
        #   • and the core should contain exactly one oxygen.
        if (thf_count == 1 and total_rings in (3, 4) and 
            (carbocyclic_count == total_rings - 1) and
            14 <= core_carbons <= 16 and core_oxygens == 1):
            reason = ("Fused ring system with {} rings ({} THF and {} carbocyclic) detected having {} core carbons and {} core oxygen, "
                      "consistent with a dihydroagarofuran sesquiterpenoid skeleton."
                      .format(total_rings, thf_count, carbocyclic_count, core_carbons, core_oxygens))
            return True, reason

    return False, "No fused ring system consistent with a dihydroagarofuran skeleton was found"

# For testing when run as main
if __name__ == "__main__":
    # Example: wilforidine (a known dihydroagarofuran sesquiterpenoid)
    test_smiles = "[C@@]12([C@@H]([C@@H]([C@@]3([C@H]([C@]14[C@]([C@H]([C@@H]([C@@H]2OC(=O)C)O)OC(C(CCC5=NC=CC=C5C(OC[C@@]3(O4)C)=O)(C)O)=O)(O)C)OC(C)=O)[H])OC(C)=O)OC(=O)C)COC(C)=O"
    result, explanation = is_dihydroagarofuran_sesquiterpenoid(test_smiles)
    print("Result:", result)
    print("Explanation:", explanation)