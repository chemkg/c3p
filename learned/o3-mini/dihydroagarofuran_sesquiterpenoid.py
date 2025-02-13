"""
Classifies: CHEBI:71548 dihydroagarofuran sesquiterpenoid
"""
#!/usr/bin/env python
"""
Classifies: dihydroagarofuran sesquiterpenoid
Definition: Any sesquiterpenoid with a dihydroagarofuran skeleton.

Heuristic:
  1. Look for at least one tetrahydrofuran ring (5-membered ring with 1 oxygen and 4 carbons).
  2. Check if that ring is fused with at least two other rings (i.e. shares at least 2 atoms with each additional ring).
  3. In the fused ring system, roughly 10–20 carbons should be present (the true dihydroagarofuran core has 15 carbons).
If these conditions are met then we classify the molecule as a dihydroagarofuran sesquiterpenoid.
"""
from rdkit import Chem

def is_dihydroagarofuran_sesquiterpenoid(smiles: str):
    """
    Determines if a molecule is a dihydroagarofuran sesquiterpenoid based on its SMILES string.
    
    This function uses a heuristic:
      - It searches for a tetrahydrofuran ring (a 5-membered ring with one oxygen and four carbons).
      - It then checks that this ring is fused (sharing at least two atoms) with at least two other rings.
      - It also checks that the union of atoms in the fused ring system has a carbon count roughly 
        consistent with a sesquiterpenoid core (about 15 carbons, here allowed in the 10 to 20 range).
        
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a dihydroagarofuran sesquiterpenoid, False otherwise.
        str: Explanation for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get all rings in the molecule as tuples of atom indices.
    ring_info = mol.GetRingInfo()
    if not ring_info or not ring_info.AtomRings():
        return False, "No rings detected in the molecule"
    rings = ring_info.AtomRings()
    
    # Look for candidate tetrahydrofuran rings:
    candidate_thf_rings = []
    for ring in rings:
        if len(ring) == 5:
            oxygen_count = 0
            carbon_count = 0
            for idx in ring:
                symbol = mol.GetAtomWithIdx(idx).GetSymbol()
                if symbol == 'O':
                    oxygen_count += 1
                elif symbol == 'C':
                    carbon_count += 1
            # A tetrahydrofuran ring typically has 1 O and 4 C.
            if oxygen_count == 1 and carbon_count == 4:
                candidate_thf_rings.append(set(ring))
    
    if not candidate_thf_rings:
        return False, "No tetrahydrofuran (5-membered, 1O/4C) ring found"
    
    # Now, for each candidate THF ring, check for fusion with other rings:
    for thf_ring in candidate_thf_rings:
        fused_rings = [thf_ring]  # start with the THF ring
        # Look for rings that share at least 2 atoms with our candidate THF ring.
        for other_ring in rings:
            other_ring_set = set(other_ring)
            if other_ring_set == thf_ring:
                continue
            if len(thf_ring.intersection(other_ring_set)) >= 2:
                fused_rings.append(other_ring_set)
        if len(fused_rings) < 3:
            # Not enough fused rings to match the dihydroagarofuran skeleton.
            continue
        
        # Compute union of atoms in the fused ring system.
        union_atoms = set()
        for ring_set in fused_rings:
            union_atoms = union_atoms.union(ring_set)
        # Count carbons in the fused system.
        core_carbons = 0
        for idx in union_atoms:
            if mol.GetAtomWithIdx(idx).GetSymbol() == 'C':
                core_carbons += 1
        # Many dihydroagarofuran cores are sesquiterpenoids (roughly 15 carbons).
        # Allow some variation with substituents removed (here we require between 10 and 20 C).
        if 10 <= core_carbons <= 20:
            return True, ("Fused ring system detected with a tetrahydrofuran ring and ~{} core carbons (expected "
                          "for a sesquiterpenoid).".format(core_carbons))
    
    return False, "No fused ring system consistent with a dihydroagarofuran skeleton was found"

# If run as a script, some test cases can be executed.
if __name__ == "__main__":
    # Example: one of the provided SMILES for wilforidine.
    test_smiles = "[C@@]12([C@@H]([C@@H]([C@@]3([C@H]([C@]14[C@]([C@H]([C@@H]([C@@H]2OC(=O)C)O)OC(C(CCC5=NC=CC=C5C(OC[C@@]3(O4)C)=O)(C)O)=O)(O)C)OC(C)=O)[H])OC(C)=O)OC(=O)C)COC(C)=O"
    result, reason = is_dihydroagarofuran_sesquiterpenoid(test_smiles)
    print("Result:", result)
    print("Reason:", reason)