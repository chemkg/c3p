"""
Classifies: CHEBI:48953 cyclohexenones
"""
"""
Classifies: Cyclohexenones – any six-membered alicyclic ketone having one double bond in the ring.
The molecule must contain a six‐membered, non‐aromatic (carbon only) ring that has exactly one internal double bond,
and exactly one ketone group defined as an exocyclic double bond O attached to a ring carbon. The ketone must be conjugated,
i.e. its ring carbon shares a bond with one of the two atoms forming the ring double bond.
"""

from rdkit import Chem

def is_cyclohexenones(smiles: str):
    """
    Determines if a molecule is a cyclohexenone based on its SMILES string.
    
    A cyclohexenone is defined as a six‐membered (non‐aromatic) carbon ring that contains exactly one internal double bond,
    and exactly one ketone functional group (a ring carbon double-bonded to an exocyclic oxygen) such that
    the ketone carbon is adjacent to one of the ring carbons forming the double bond.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule meets the cyclohexenone criteria, False otherwise.
        str: Reason for classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()  # each ring is a tuple of atom indices
    
    # Loop over candidate rings
    for ring in rings:
        if len(ring) != 6:
            continue  # only interested in six-membered rings
        
        # Check that every atom in the ring is carbon and is non-aromatic.
        ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
        if any(atom.GetAtomicNum() != 6 for atom in ring_atoms):
            continue
        if any(atom.GetIsAromatic() for atom in ring_atoms):
            continue
            
        ring_set = set(ring)
        
        # Identify double bonds within the ring.
        double_bonds = []
        for bond in mol.GetBonds():
            # Both atoms in the bond must belong to the ring.
            a1 = bond.GetBeginAtomIdx()
            a2 = bond.GetEndAtomIdx()
            if a1 in ring_set and a2 in ring_set:
                if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                    # To avoid double-counting, store sorted tuple.
                    bond_tuple = tuple(sorted([a1,a2]))
                    if bond_tuple not in double_bonds:
                        double_bonds.append(bond_tuple)
        
        if len(double_bonds) != 1:
            # We require exactly one double bond in the ring.
            continue
        
        # Search for ketone groups in the ring.
        # A ketone group here means one ring carbon has a double bond to an oxygen (external to the ring).
        ketone_candidates = []
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            # Examine neighbors that are not in the ring.
            for nbr in atom.GetNeighbors():
                if nbr.GetAtomicNum() == 8 and nbr.GetIdx() not in ring_set:
                    bond = mol.GetBondBetweenAtoms(idx, nbr.GetIdx())
                    if bond and bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                        ketone_candidates.append(idx)
                        
        # We require exactly one ketone group attached to the ring.
        if len(ketone_candidates) != 1:
            continue
        
        ketone_atom_idx = ketone_candidates[0]
        # Enforce conjugation: the ketone carbon in the ring must be adjacent (within the ring) to one of the double-bond atoms.
        dbond_a, dbond_b = double_bonds[0]
        # Get ring neighbors (by bond connectivity) of the ketone atom within the ring.
        ketone_neighbors = [nbr.GetIdx() for nbr in mol.GetAtomWithIdx(ketone_atom_idx).GetNeighbors() if nbr.GetIdx() in ring_set]
        if dbond_a in ketone_neighbors or dbond_b in ketone_neighbors:
            return True, "Found a cyclohexenone: six-membered non-aromatic carbon ring with one double bond and conjugated ketone"
        else:
            # The ketone group is not conjugated with the double bond.
            continue
            
    return False, "No cyclohexenone ring (six-membered, non-aromatic, with one double bond and conjugated ketone) found"

# Example usage: testing with some candidate SMILES, including the simple cyclohex-2-enone.
if __name__ == "__main__":
    test_smiles = [
        "O=C1CCCC=C1",         # cyclohex-2-enone, should be True
        "O=C1CC=CCC1",         # alternative depiction, should be True
        "O=C1C=CCCC1"          # slight connectivity variation, should be True if it meets criteria
    ]
    for smi in test_smiles:
        result, reason = is_cyclohexenones(smi)
        print(f"SMILES: {smi}\nClassified as cyclohexenone: {result}\nReason: {reason}\n")