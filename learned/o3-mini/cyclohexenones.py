"""
Classifies: CHEBI:48953 cyclohexenones
"""
"""
Classifies: Cyclohexenones – Any six-membered alicyclic ketone having one double bond in the ring.
Improved version: In addition to ring size, non-aromaticity, exactly one ring double bond, and a C=O (ketone) group, 
we require that the ketone carbon is conjugated with the ring double bond (i.e. it is directly bonded to one of the double‐bonded ring carbons).
"""

from rdkit import Chem

def is_cyclohexenones(smiles: str):
    """
    Determines if a molecule is a cyclohexenone based on its SMILES string.
    A cyclohexenone is defined as a six-membered alicyclic (non-aromatic) carbon ring that contains exactly one internal double bond,
    and exactly one ketone functional group (a ring carbon double-bonded to an exocyclic oxygen) such that the carbonyl is conjugated 
    with the double bond (i.e. the ketone carbon is adjacent in the ring to one of the carbons involved in the ring double bond).
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule meets the cyclohexenone criteria, False otherwise.
        str: Reason for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get ring information from the molecule.
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()  # each ring is a tuple of atom indices
    
    # Iterate over rings looking for a candidate six-membered ring.
    for ring in rings:
        if len(ring) != 6:
            continue  # only interested in six-membered rings
            
        # Check that every atom in the ring is a carbon and non-aromatic.
        ring_atoms = [mol.GetAtomWithIdx(i) for i in ring]
        if any(atom.GetAtomicNum() != 6 for atom in ring_atoms):
            continue
        if any(atom.GetIsAromatic() for atom in ring_atoms):
            continue
        
        # Prepare a lookup for quick membership tests.
        ring_set = set(ring)
        
        # Count double bonds within the ring and record which bond indices (atom pairs) are double.
        double_bonds = []
        for i in ring:
            atom = mol.GetAtomWithIdx(i)
            for bond in atom.GetBonds():
                a1 = bond.GetBeginAtomIdx()
                a2 = bond.GetEndAtomIdx()
                # Consider bonds with both atoms in the ring (avoid counting an external bond)
                if a1 in ring_set and a2 in ring_set:
                    # To avoid double-counting, ensure we only consider bond once.
                    if a1 < a2 and bond.GetBondTypeAsDouble() == 2:
                        double_bonds.append((a1, a2))
        
        if len(double_bonds) != 1:
            continue  # must have exactly one double bond in the ring
        
        # Look for exactly one ketone group in the ring.
        # A ketone here is defined as a ring carbon that is double-bonded to an oxygen (the O must be outside the ring).
        ketone_atom = None
        ketone_count = 0
        for i in ring:
            atom = mol.GetAtomWithIdx(i)
            # Examine neighbors, but only consider exocyclic oxygens.
            for nbr in atom.GetNeighbors():
                if nbr.GetAtomicNum() == 8 and nbr.GetIdx() not in ring_set:
                    bond = mol.GetBondBetweenAtoms(i, nbr.GetIdx())
                    if bond and bond.GetBondTypeAsDouble() == 2:
                        ketone_atom = i
                        ketone_count += 1
            # Early exit if more than one ketone is found in the ring.
            if ketone_count > 1:
                break
        if ketone_count != 1:
            continue  # candidate ring must have exactly one ketone group
        
        # Now enforce conjugation:
        # The ketone carbon must be adjacent (in the ring) to one of the atoms that form the ring double bond.
        # First, unpack the double bond pair.
        dbond_a, dbond_b = double_bonds[0]
        # Get neighbors of the ketone atom that are in the ring.
        ketone_neighbors = [nbr.GetIdx() for nbr in mol.GetAtomWithIdx(ketone_atom).GetNeighbors() if nbr.GetIdx() in ring_set]
        if dbond_a in ketone_neighbors or dbond_b in ketone_neighbors:
            return True, "Found a six-membered non-aromatic ring with one double bond and a ketone group conjugated to it"
        else:
            # The ring had the requisite groups but the ketone is not conjugated with the double bond.
            continue
            
    return False, "No cyclohexenone ring (six-membered, non-aromatic, with one double bond and a conjugated ketone) found"

# Example usage:
if __name__ == "__main__":
    # Test with a few SMILES strings (the simplest being cyclohex-2-enone)
    test_smiles = [
        "O=C1CCCC=C1",  # cyclohex-2-enone (should be True)
        "O=C1CC=CCC1",  # alternative depiction of cyclohex-2-enone (should be True)
        "O=C1C=CCCC1"   # structure with same formulas but maybe different connectivity
    ]
    for smi in test_smiles:
        result, reason = is_cyclohexenones(smi)
        print(f"SMILES: {smi}\nClassified as cyclohexenone: {result}\nReason: {reason}\n")