"""
Classifies: CHEBI:71548 dihydroagarofuran sesquiterpenoid
"""
#!/usr/bin/env python
"""
Classifies: dihydroagarofuran sesquiterpenoid
Definition: Any sesquiterpenoid with a dihydroagarofuran skeleton.
Heuristic (improved):
  1. Search for candidate tetrahydrofuran (THF) rings:
     A THF candidate is defined as a 5-membered ring containing exactly one oxygen atom and 4 carbons.
  2. For each candidate THF ring, build the fused ring system (i.e. include any ring sharing at least one atom).
  3. In that fused system, count:
       - the number of rings that exactly match the THF criteria (should be exactly one)
       - the number of rings that are fully carbocyclic (all atoms are carbon).
  4. Also count the total number of carbons in the union of all atoms of the fused rings.
  5. If the fused system contains exactly one THF and two to three carbocyclic rings (i.e. usually a 3–4 ring system),
     and the fused system contains 14–16 carbons (consistent with a 15-carbon sesquiterpenoid core),
     then we classify the molecule as a dihydroagarofuran sesquiterpenoid.
    
Note: This heuristic is approximate and may not correctly classify every structure.
"""

from rdkit import Chem

def is_dihydroagarofuran_sesquiterpenoid(smiles: str):
    """
    Determines if a molecule is a dihydroagarofuran sesquiterpenoid based on its SMILES string.
    
    The function uses an improved heuristic:
      - It looks for candidate tetrahydrofuran rings (5-membered rings with 1 oxygen and 4 carbons).
      - For each candidate, it gathers all rings fused (sharing at least one atom) with the candidate.
      - It then requires that in that fused system:
            * Exactly one ring meets the strict THF criteria.
            * The remainder of rings are entirely carbocyclic.
            * The total number of carbon atoms in the union of the fused rings is in the range 14–16.
      - If these conditions are met, the molecule is classified as a dihydroagarofuran sesquiterpenoid.
    
    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if molecule is classified as a dihydroagarofuran sesquiterpenoid, False otherwise.
        str: Explanation for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    ring_info = mol.GetRingInfo()
    if not ring_info or not ring_info.AtomRings():
        return False, "No rings detected in the molecule"
    rings = ring_info.AtomRings()  # each ring as a tuple of atom indices
    
    # Helper function: check if a ring is a candidate THF.
    def is_THF(ring):
        if len(ring) != 5:
            return False
        oxygen_count = 0
        carbon_count = 0
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            sym = atom.GetSymbol()
            if sym == "O":
                oxygen_count += 1
            elif sym == "C":
                carbon_count += 1
            else:
                # any extra heteroatom disqualifies THF for our purpose
                return False
        return (oxygen_count == 1 and carbon_count == 4)
    
    # Identify candidate THF rings
    candidate_thf = []
    for ring in rings:
        if is_THF(ring):
            candidate_thf.append(set(ring))
    
    if not candidate_thf:
        return False, "No tetrahydrofuran (5-membered, 1O/4C) ring found"
    
    # Helper to build fused ring system starting from one candidate ring.
    # Two rings are fused if they share at least one atom.
    def fused_ring_system(start_ring):
        fused = [start_ring]
        changed = True
        while changed:
            changed = False
            for ring in rings:
                ring_set = set(ring)
                # if already in fused (based on set-equality), skip
                if any(ring_set == r for r in fused):
                    continue
                if any(len(ring_set.intersection(r)) >= 1 for r in fused):
                    fused.append(ring_set)
                    changed = True
        return fused

    # Analyze each candidate THF ring’s fused system:
    for candidate in candidate_thf:
        fused_rings = fused_ring_system(candidate)
        # Count rings by type in the fused system:
        thf_count = 0
        carbocyclic_count = 0
        for ring_set in fused_rings:
            # If ring qualifies as THF, count as such.
            if is_THF(list(ring_set)):
                thf_count += 1
            else:
                # Check if every atom in the ring is carbon
                if all(mol.GetAtomWithIdx(idx).GetSymbol() == "C" for idx in ring_set):
                    carbocyclic_count += 1
        # Build the union of all atoms in the fused system:
        fused_atom_union = set()
        for ring_set in fused_rings:
            fused_atom_union |= ring_set
        core_carbons = sum(1 for idx in fused_atom_union if mol.GetAtomWithIdx(idx).GetSymbol() == "C")
        
        # We expect:
        #   • exactly one THF ring,
        #   • and typically 2 (or sometimes 3) fully carbocyclic rings (so fused system size of 3–4 rings),
        #   • and the total carbon count of the fused system in the range 14–16.
        if (thf_count == 1 and (carbocyclic_count == 2 or carbocyclic_count == 3)
            and 14 <= core_carbons <= 16):
            reason = ("Fused ring system with {} rings ({} THF and {} carbocyclic) detected having {} core carbons, "
                      "consistent with a dihydroagarofuran sesquiterpenoid skeleton."
                      .format(len(fused_rings), thf_count, carbocyclic_count, core_carbons))
            return True, reason
    
    return False, "No fused ring system consistent with a dihydroagarofuran skeleton was found"

# When run as a script, an example test with a known dihydroagarofuran (wilforidine) is executed.
if __name__ == "__main__":
    # Example SMILES for wilforidine (a known dihydroagarofuran sesquiterpenoid):
    test_smiles = "[C@@]12([C@@H]([C@@H]([C@@]3([C@H]([C@]14[C@]([C@H]([C@@H]([C@@H]2OC(=O)C)O)OC(C(CCC5=NC=CC=C5C(OC[C@@]3(O4)C)=O)(C)O)=O)(O)C)OC(C)=O)[H])OC(C)=O)OC(=O)C)COC(C)=O"
    result, explanation = is_dihydroagarofuran_sesquiterpenoid(test_smiles)
    print("Result:", result)
    print("Explanation:", explanation)