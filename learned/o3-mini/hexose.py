"""
Classifies: CHEBI:18133 hexose
"""
"""
Classifies: Hexose
Defined as: Any six‐carbon monosaccharide which in its linear form contains either 
an aldehyde group at position 1 (aldohexose) or a ketone group at position 2 (ketohexose).
This implementation first checks for open‐chain hexose motifs using SMARTS.
If that fails, it then looks for cyclic sugar rings (either pyranose or furanose forms)
by scanning the ring systems, checking for the expected number of oxygens inside the ring
and for exocyclic CH2OH groups that complete a six‐carbon unit.
"""

from rdkit import Chem

def is_hexose(smiles: str):
    """
    Determines if a molecule contains a hexose moiety (six‐carbon monosaccharide)
    based on its SMILES string. It attempts two strategies:
      (1) Look for an open‐chain aldohexose or ketohexose SMARTS pattern.
      (2) Look for cyclic sugar rings having the expected ring size and one oxygen
          plus exocyclic –CH2OH group(s) so that six carbons are contained.
          
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if a hexose-like moiety is found, False otherwise.
        str: Explanation for classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # ---- Strategy 1: Open-chain hexose patterns ----
    # We define SMARTS patterns (ignoring stereochemistry) for a typical aldohexose open chain.
    # The pattern here requires six carbons in sequence with hydroxyl groups on C2–C6.
    # (For an aldohexose the first carbon is a formyl group.)
    aldo_hex_smarts = "[CX3H1](=O)[CX4H](O)[CX4H](O)[CX4H](O)[CX4H](O)[CH2]O"
    keto_hex_smarts = "O[C][CX4H](O)[CX4H](O)[CX4H](O)[CX4H](O)C(=O)C"
    # Note: the ketose pattern is defined with an -O at the beginning (for the CH2OH group)
    # so that the carbonyl ends up at position 2.
    try:
        aldo_mol = Chem.MolFromSmarts(aldo_hex_smarts)
        keto_mol = Chem.MolFromSmarts(keto_hex_smarts)
    except Exception as e:
        return False, f"Error in SMARTS parsing: {str(e)}"
    
    if mol.HasSubstructMatch(aldo_mol):
        return True, "Aldehyde-based open-chain hexose motif detected; consistent with an aldohexose"
    if mol.HasSubstructMatch(keto_mol):
        return True, "Ketone-based open-chain hexose motif detected; consistent with a ketohexose"
    
    # ---- Strategy 2: Cyclic (ring) hexose detection ----
    # For cyclic sugars the hexose motif usually appears as a ring with one heteroatom (O)
    # plus (for pyranoses) one exocyclic CH2OH branch (5 ring carbons + 1 extra carbon = 6 carbons)
    # or (for furanoses) a 5-membered ring (4 ring carbons + 2 extra carbons = 6 carbons).
    ring_info = mol.GetRingInfo()
    if not ring_info:
        return False, "No rings found and no open-chain hexose motif detected"
    
    for ring in ring_info.AtomRings():
        ring_size = len(ring)
        # Get ring atoms
        ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
        # Count oxygen atoms in the ring
        oxygens_in_ring = sum(1 for atom in ring_atoms if atom.GetAtomicNum() == 8)
        
        # Pyranose candidate: 6-membered (5 ring carbons, one oxygen)
        if ring_size == 6 and oxygens_in_ring == 1:
            # We now try to see if exactly one exocyclic carbon (CH2OH) is attached to the ring,
            # so that total C count = (5 in ring + 1 exocyclic) = 6.
            exo_c_count = 0
            for idx in ring:
                atom = mol.GetAtomWithIdx(idx)
                if atom.GetAtomicNum() == 6:
                    # Look at neighbors that are not in the ring.
                    for nb in atom.GetNeighbors():
                        if nb.GetIdx() not in ring and nb.GetAtomicNum() == 6:
                            # Check if this neighbor is a CH2 group: exactly 2 hydrogens and a single bond to an O
                            # (the oxygen part for –CH2OH)
                            if nb.GetTotalNumHs() >= 2:
                                # Check if one of its neighbors (not in the ring) is an -OH:
                                for nb2 in nb.GetNeighbors():
                                    if nb2.GetIdx() not in ring and nb2.GetAtomicNum() == 8:
                                        # Found a candidate exocyclic CH2OH.
                                        exo_c_count += 1
                                        break
            if exo_c_count >= 1:
                return True, f"Cyclic (pyranose) sugar motif detected (ring size 6 with one oxygen and exocyclic CH2OH group)"
        
        # Furanose candidate: 5-membered ring expected (4 carbons + 1 oxygen).
        # In many hexoses that adopt a furanose form there are 2 exocyclic carbons (to reach total C=6).
        if ring_size == 5 and oxygens_in_ring == 1:
            exo_c_count = 0
            for idx in ring:
                atom = mol.GetAtomWithIdx(idx)
                if atom.GetAtomicNum() == 6:
                    for nb in atom.GetNeighbors():
                        if nb.GetIdx() not in ring and nb.GetAtomicNum() == 6:
                            # Very roughly check for CH2OH (at least 1 or 2 such branches).
                            if nb.GetTotalNumHs() >= 2:
                                for nb2 in nb.GetNeighbors():
                                    if nb2.GetIdx() not in ring and nb2.GetAtomicNum() == 8:
                                        exo_c_count += 1
                                        break
            # For furanose form, we expect about 2 exocyclic carbons attached.
            if exo_c_count >= 2:
                return True, f"Cyclic (furanose) sugar motif detected (ring size 5 with one oxygen and two exocyclic CH2OH-like groups)"
    
    return False, "No hexose substructure detected (neither open-chain carbonyl patterns nor cyclic sugar motif found)"

# Example usage (you can remove or comment out these lines when incorporating the function into another program)
if __name__ == "__main__":
    # Try one of the provided examples: alpha-D-glucose
    test_smiles = "OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@@H]1O"
    result, reason = is_hexose(test_smiles)
    print("Test result:", result)
    print("Reason:", reason)