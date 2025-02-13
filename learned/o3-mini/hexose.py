"""
Classifies: CHEBI:18133 hexose
"""
#!/usr/bin/env python
"""
Classifies: Hexose
Defined as: Any six‐carbon monosaccharide which in its linear form contains either 
an aldehyde group at position 1 (aldohexose) or a ketone group at position 2 (ketohexose).

This implementation uses two strategies:
  Strategy 1. Try to match open‐chain hexose SMARTS patterns for the aldo- and keto- forms.
  Strategy 2. Look for cyclic sugar substructures (pyranose or furanose) by scanning the molecule’s rings.
      For a pyranose candidate a 6-membered ring should contain 5 carbons and 1 oxygen,
      and exactly one exocyclic carbon branch (of “CH2OH” type) should be attached.
      For a furanose candidate a 5-membered ring with 4 carbons and 1 oxygen
      should have exactly two such exocyclic carbon branches.
      
If a candidate hexose substructure is found, True is returned along with a reason;
otherwise, False is returned.
"""

from rdkit import Chem

def is_hexose(smiles: str):
    """
    Determines if a molecule contains a hexose moiety (six‐carbon monosaccharide)
    based on its SMILES string. Two strategies are attempted:
      1. Open-chain detection: using SMARTS patterns for aldohexose and ketohexose.
      2. Cyclic detection: scanning each ring to see if carbon count (ring carbons plus 
         exocyclic CH2OH-like substitutions) equals six.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if a hexose-like moiety is found, False otherwise.
        str: Explanation for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # -------- Strategy 1: Open-chain detection --------
    # We require exactly six carbons with the appropriate terminal carbonyl.
    # Aldohexose: must start with a formyl group and be followed by four CH-O groups and finish with CH2OH.
    # (SMARTS pattern written without stereochemistry for generality)
    aldo_hex_smarts = "[CX3H1](=O)[CX4](O)[CX4](O)[CX4](O)[CX4](O)[CH2]O"
    # Ketohexose: HOCH2- at the beginning, then a CH-O group, then a carbonyl at C3 (which will be position 2 in the chain),
    # followed by two CH-O groups and finishing with CH2OH.
    keto_hex_smarts = "[CH2]O[CX4](O)[CX4](=O)[CX4](O)[CX4](O)[CH2]O"
    
    try:
        aldo_pattern = Chem.MolFromSmarts(aldo_hex_smarts)
        keto_pattern = Chem.MolFromSmarts(keto_hex_smarts)
    except Exception as e:
        return False, f"Error parsing SMARTS: {str(e)}"
    
    # Check open-chain aldohexose match.
    matches = mol.GetSubstructMatches(aldo_pattern)
    if matches:
        # Optional: you might check the matching tuple to ensure six carbons, etc.
        return True, "Aldehyde-based open-chain hexose motif detected; consistent with an aldohexose"

    # Check open-chain ketohexose match.
    matches = mol.GetSubstructMatches(keto_pattern)
    if matches:
        return True, "Ketone-based open-chain hexose motif detected; consistent with a ketohexose"
    
    # -------- Strategy 2: Cyclic detection --------
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    if not rings:
        return False, "No rings found and no open-chain hexose motif detected"
    
    # Helper function to test if an atom is a potential CH2OH fragment:
    def is_ch2oh(candidate):
        # Must be a carbon with exactly two hydrogens (we allow >=2 since implicit H's might be counted)
        # and at least one neighbor (outside the candidate substructure) is oxygen.
        if candidate.GetAtomicNum() != 6:
            return False
        if candidate.GetTotalNumHs() < 2:
            return False
        for nb in candidate.GetNeighbors():
            if nb.GetAtomicNum() == 8:
                return True
        return False

    # Iterate over each ring:
    for ring in rings:
        ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
        ring_size = len(ring)
        # Count number of oxygens and carbons within the ring.
        oxygens_in_ring = sum(1 for atom in ring_atoms if atom.GetAtomicNum() == 8)
        carbons_in_ring = sum(1 for atom in ring_atoms if atom.GetAtomicNum() == 6)
        
        # ----- Check for pyranose candidate -----
        # A typical pyranose ring: 6 atoms with exactly one oxygen and 5 carbons.
        if ring_size == 6 and oxygens_in_ring == 1 and carbons_in_ring == 5:
            # For a hexose, we expect one exocyclic carbon (CH2OH) attached to one of the ring carbons.
            exo_carbon_ids = set()
            for idx in ring:
                atom = mol.GetAtomWithIdx(idx)
                # Only consider carbon atoms in the ring.
                if atom.GetAtomicNum() != 6:
                    continue
                # Look for neighbors outside the ring.
                for nb in atom.GetNeighbors():
                    if nb.GetIdx() not in ring and nb.GetAtomicNum() == 6:
                        if is_ch2oh(nb):
                            exo_carbon_ids.add(nb.GetIdx())
            if len(exo_carbon_ids) == 1:
                return True, "Cyclic (pyranose) sugar motif detected: 6-membered ring (5 carbons + 1 O) with one exocyclic CH2OH branch forming 6 carbons"
        
        # ----- Check for furanose candidate -----
        # A typical furanose ring: 5 atoms with exactly one oxygen and 4 carbons.
        if ring_size == 5 and oxygens_in_ring == 1 and carbons_in_ring == 4:
            # For a hexose in a furanose form, we expect two exocyclic carbons (each as CH2OH) so that 4 ring carbons + 2 exo = 6.
            exo_carbon_ids = set()
            for idx in ring:
                atom = mol.GetAtomWithIdx(idx)
                if atom.GetAtomicNum() != 6:
                    continue
                for nb in atom.GetNeighbors():
                    if nb.GetIdx() not in ring and nb.GetAtomicNum() == 6:
                        if is_ch2oh(nb):
                            exo_carbon_ids.add(nb.GetIdx())
            if len(exo_carbon_ids) == 2:
                return True, "Cyclic (furanose) sugar motif detected: 5-membered ring (4 carbons + 1 O) with two exocyclic CH2OH branches forming 6 carbons"
    
    return False, "No hexose substructure detected (neither open-chain patterns nor cyclic sugar motif matched)"

# Example usage:
if __name__ == "__main__":
    # You can test with one of the provided examples such as alpha-D-glucose.
    test_smiles = "OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@@H]1O"
    result, reason = is_hexose(test_smiles)
    print("Test result:", result)
    print("Reason:", reason)