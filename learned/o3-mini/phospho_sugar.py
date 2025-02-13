"""
Classifies: CHEBI:33447 phospho sugar
"""
"""
Classifies: Phospho sugars 
Definition: Any monosaccharide (here assumed to contain a 5- or 6-membered ring with exactly one ring oxygen) 
that contains an alcoholic –OH that is esterified with phosphoric acid (i.e. has a phosphate group attached via oxygen).
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_phospho_sugar(smiles: str):
    """
    Determines if a molecule is a phospho sugar based on its SMILES string.
    Definition (in this implementation):
      - The molecule contains a sugar ring, defined as a 5- or 6-membered ring with exactly one oxygen atom.
      - The molecule contains a phosphate ester group. This is detected by finding a phosphorus atom
        that has at least one double-bonded oxygen (indicating a P=O) and that is connected via an oxygen 
        (the ester oxygen) to either the ring directly or to a substituent attached to a ring atom.
    
    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if molecule is classified as a phospho sugar, False otherwise.
        str: Reason for the classification decision.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # --- STEP 1: Find a candidate sugar ring ---
    # We look for rings of size 5 or 6 that have exactly one oxygen atom.
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    sugar_ring_indices = []  # will hold indices of atoms in one of the sugar rings, if found
    found_sugar_ring = False
    
    for ring in atom_rings:
        if len(ring) in (5, 6):
            # Count the number of oxygen atoms in the ring.
            o_count = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8)
            if o_count == 1:
                sugar_ring_indices = list(ring)
                found_sugar_ring = True
                break
    if not found_sugar_ring:
        return False, "No sugar ring detected (no 5- or 6-membered ring with exactly one oxygen found)"

    # --- STEP 2: Look for a phosphate ester group ---
    # A phosphate group typically has a phosphorus atom (atomic num 15)
    # bonded to at least one oxygen atom via a double bond (P=O).
    phosphate_found = False
    # Also track if the connecting oxygen is linked to the sugar ring.
    phosphate_connected_to_sugar = False

    # Iterate over all atoms looking for an oxygen that is connected to a phosphorus.
    for atom in mol.GetAtoms():
        # Check if current atom is oxygen (atomic num 8)
        if atom.GetAtomicNum() != 8:
            continue
        for neighbor in atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 15:  # phosphorus found
                p_atom = neighbor
                # Check that the phosphorus has at least one double-bonded oxygen (P=O)
                has_double_bonded_O = False
                for bond in p_atom.GetBonds():
                    # Get the other atom of the bond.
                    other = bond.GetOtherAtom(p_atom)
                    if other.GetAtomicNum() == 8 and bond.GetBondTypeAsDouble() == 2.0:
                        has_double_bonded_O = True
                        break
                if not has_double_bonded_O:
                    continue  # not a phosphate-like group

                phosphate_found = True
                # Now, we expect the O (atom) that is attached to phosphorus to be an ester oxygen.
                # Check if this oxygen is either part of the sugar ring or attached to an atom in the sugar ring.
                # First, check if the oxygen itself is in sugar ring.
                if atom.GetIdx() in sugar_ring_indices:
                    phosphate_connected_to_sugar = True
                    break
                # Otherwise, check its neighbors: is one of them in the sugar ring?
                for o_neighbor in atom.GetNeighbors():
                    if o_neighbor.GetIdx() in sugar_ring_indices:
                        phosphate_connected_to_sugar = True
                        break
                if phosphate_connected_to_sugar:
                    break
        if phosphate_found and phosphate_connected_to_sugar:
            break

    if not phosphate_found:
        return False, "No phosphate group (with a P=O bond) found in the molecule"
    
    if not phosphate_connected_to_sugar:
        return False, "Phosphate group found but not connected to the detected sugar ring"

    # --- If both criteria are met, classify as phospho sugar ---
    return True, "Molecule contains a sugar ring and a phosphate ester group attached via an alcoholic oxygen"

# (Optional) For testing, one can run:
if __name__ == "__main__":
    test_smiles = "O[C@H](COP(O)(O)=O)[C@@H](O)CC=O"  # 2-deoxy-D-ribose 5-phosphate example
    result, reason = is_phospho_sugar(test_smiles)
    print("Result:", result)
    print("Reason:", reason)