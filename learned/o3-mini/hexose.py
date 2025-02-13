"""
Classifies: CHEBI:18133 hexose
"""
#!/usr/bin/env python
"""
Classifies: Hexose
Defined as: Any six‐carbon monosaccharide which in its linear form contains either an aldehyde group at position 1 (aldohexose)
or a ketone group at position 2 (ketohexose). In the cyclic form, a hexose is expected to have a ring that, when combined with
its exocyclic CH2OH substituents, accounts for exactly six carbons.
  
This improved implementation uses two strategies:
  Strategy 1. Open‐chain detection via SMARTS patterns (for aldo- and ketohexoses).
  Strategy 2. Cyclic detection: examine each ring that contains exactly one oxygen and count the carbons in the ring plus
              the number of exocyclic CH2OH-like substituents. For a pyranose, a six-membered ring should have 5 ring carbons
              and one CH2OH branch; for a furanose, a five-membered ring (4 carbons) should have two such branches.
If a candidate hexose substructure is found, the function returns True and explains its decision.

Note: When the hexose motif is part of a larger molecule (e.g. in glycosides), the algorithm may flag
it as hexose because it looks like a free hexose moiety.
"""

from rdkit import Chem

def is_hexose(smiles: str):
    """
    Determines if a molecule contains a hexose moiety (six‐carbon monosaccharide) based on its SMILES string.
    Two strategies are used:
      1. Open-chain detection: using SMARTS patterns for aldohexose and ketohexose.
      2. Cyclic detection: scanning rings for typical sugar rings (pyranose or furanose) by examining both
         atoms in the ring and out-of-ring (exocyclic) CH2OH substituents.
         
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule contains a hexose-like motif, False otherwise.
        str: Explanation of the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # -------- Strategy 1: Open-chain detection --------
    # Here we try to match open-chain aldohexose and ketohexose patterns.
    # Aldohexose: expected to start with an aldehyde then four CH-O units and finish with CH2OH.
    aldo_hex_smarts = "[CX3H1](=O)[CX4](O)[CX4](O)[CX4](O)[CX4](O)[CH2]O"
    # Ketohexose: might have HOCH2- then CH-O, then a carbonyl (ketone) then two CH-O groups and CH2OH.
    keto_hex_smarts = "[CH2]O[CX4](O)[CX4](=O)[CX4](O)[CX4](O)[CH2]O"

    try:
        aldo_pattern = Chem.MolFromSmarts(aldo_hex_smarts)
        keto_pattern = Chem.MolFromSmarts(keto_hex_smarts)
    except Exception as e:
        return False, f"Error parsing SMARTS: {str(e)}"
    
    if mol.HasSubstructMatch(aldo_pattern):
        return True, "Aldehyde-based open-chain hexose motif detected; consistent with an aldohexose"
    if mol.HasSubstructMatch(keto_pattern):
        return True, "Ketone-based open-chain hexose motif detected; consistent with a ketohexose"

    # -------- Strategy 2: Cyclic detection --------
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    if not rings:
        return False, "No rings found and no open-chain hexose motif detected"

    # Helper: Check if a candidate atom appears to be a CH2OH group.
    def is_ch2oh(candidate):
        # Candidate must be carbon and ideally should have exactly two hydrogens.
        if candidate.GetAtomicNum() != 6:
            return False
        # Using GetTotalNumHs may yield a value >=2. We assume a CH2OH carbon should have at least 2 H's.
        if candidate.GetTotalNumHs() < 2:
            return False
        # Must be attached (besides the one connection into the ring) to an oxygen atom.
        # We allow additional bonds (e.g. the OH group) but at least one neighbor (other than the ring atom) should be oxygen.
        for nb in candidate.GetNeighbors():
            if nb.GetAtomicNum() == 8:
                return True
        return False

    # Iterate over each ring to see if it could be a hexose ring.
    for ring in rings:
        # Get atoms and indices in the ring.
        ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
        # Count oxygens in the ring.
        oxygens_in_ring = sum(1 for atom in ring_atoms if atom.GetAtomicNum() == 8)
        # Proceed only if the ring has exactly one oxygen (typical for sugars).
        if oxygens_in_ring != 1:
            continue

        # For each ring, collect the indices of ring carbons.
        ring_carbons = [atom for atom in ring_atoms if atom.GetAtomicNum() == 6]

        # Based on ring size, define expected exocyclic CH2OH branches:
        # For pyranose: 6-membered ring (5 C) => expect 1 exocyclic carbon.
        # For furanose: 5-membered ring (4 C) => expect 2 exocyclic carbons.
        expected_exo = 0
        if len(ring) == 6:
            # pyranose candidate: expect exactly one exocyclic CH2OH substituent.
            expected_exo = 1
        elif len(ring) == 5:
            # furanose candidate: expect two exocyclic CH2OH substituents.
            expected_exo = 2
        else:
            continue  # Skip rings that are not typical sugar rings (not 5 or 6 members)

        exo_carbon_ids = set()
        # Check each ring carbon.
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() != 6:
                continue
            # Look for neighbors outside the ring that are carbons and look like CH2OH.
            for nb in atom.GetNeighbors():
                if nb.GetIdx() in ring:
                    continue
                if nb.GetAtomicNum() == 6 and is_ch2oh(nb):
                    exo_carbon_ids.add(nb.GetIdx())

        # Now, if the number of ring carbons plus the exocyclic carbons equals 6, we consider it hexose.
        total_carbons = len(ring_carbons) + len(exo_carbon_ids)
        if total_carbons == 6 and len(exo_carbon_ids) == expected_exo:
            if len(ring) == 6:
                return True, ("Cyclic (pyranose) sugar motif detected: 6-membered ring (5 carbons + 1 exocyclic carbon) "
                              "forming 6 total carbons")
            else:
                return True, ("Cyclic (furanose) sugar motif detected: 5-membered ring (4 carbons + 2 exocyclic carbons) "
                              "forming 6 total carbons")
    
    return False, "No hexose substructure detected (neither open-chain patterns nor cyclic sugar motif matched)"

# Example usage:
if __name__ == "__main__":
    # You can test this function with one or more provided examples.
    test_smiles = "OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@@H]1O"  # alpha-D-glucose (pyranose form)
    result, reason = is_hexose(test_smiles)
    print("Test result:", result)
    print("Reason:", reason)