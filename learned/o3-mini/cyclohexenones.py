"""
Classifies: CHEBI:48953 cyclohexenones
"""
"""
Classifies: Cyclohexenones – any six‐membered alicyclic ketone having one double bond in the ring.
A cyclohexenone is defined here as a six‐membered (non‐aromatic) ring that has exactly one internal double bond,
and at least one ring carbon bears an exocyclic ketone group (a carbon=oxygen double bond with that oxygen having no other bonds)
such that the carbonyl carbon is conjugated (i.e. directly bonded within the ring) to one of the atoms involved in the internal double bond.
We relax the requirement that every atom of the ring must be carbon in order not to miss cyclohexenone moieties in fused or substituted rings.
"""

from rdkit import Chem

def is_cyclohexenones(smiles: str):
    """
    Determines if a molecule contains a cyclohexenone moiety based on its SMILES string.
    
    The algorithm looks for a six‐membered, non‐aromatic ring (even if it contains heteroatoms) that has exactly one double bond 
    between ring atoms (an internal double bond) and at least one exocyclic ketone group attached to a ring atom. 
    In addition the carbon bearing the ketone must be adjacent (inside the ring) to one of the atoms that form the internal double bond,
    ensuring conjugation.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a cyclohexenone (i.e. it contains a cyclohexenone moiety), False otherwise.
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
            continue  # only consider six-membered rings
        
        # Require the ring to be non‐aromatic (at least no atom is marked aromatic)
        ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
        if any(atom.GetIsAromatic() for atom in ring_atoms):
            continue

        ring_set = set(ring)
        
        # Find internal double bonds in the ring (both atoms in the bond must lie in the ring)
        internal_dbonds = []
        for bond in mol.GetBonds():
            a1 = bond.GetBeginAtomIdx()
            a2 = bond.GetEndAtomIdx()
            if a1 in ring_set and a2 in ring_set:
                if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                    tup = tuple(sorted([a1, a2]))
                    if tup not in internal_dbonds:
                        internal_dbonds.append(tup)
        # We require exactly one internal double bond
        if len(internal_dbonds) != 1:
            continue

        # Identify ketone candidates on ring atoms.
        # A ketone candidate is a ring atom (preferably carbon) that has a double bond to an oxygen
        # where the oxygen is exocyclic (not in the ring) and its degree (number of bonds) is 1.
        ketone_candidates = []
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            # We require the atom to be carbon (this is typical for a cyclohexenone carbonyl carbon)
            if atom.GetAtomicNum() != 6:
                continue
            for nbr in atom.GetNeighbors():
                if nbr.GetAtomicNum() == 8 and nbr.GetIdx() not in ring_set:
                    # check that the bond between atom and nbr is a double bond
                    bond = mol.GetBondBetweenAtoms(idx, nbr.GetIdx())
                    if bond and bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                        # also check that the oxygen atom is “clean” (only bonded once)
                        if nbr.GetDegree() == 1:
                            ketone_candidates.append(idx)
        # We require at least one ketone candidate to be present.
        if not ketone_candidates:
            continue

        # For each candidate ketone, check that it is conjugated to the ring double bond.
        dbond_a, dbond_b = internal_dbonds[0]
        for ket_idx in ketone_candidates:
            # Get the ring neighbors of the ketone candidate
            ket_neighbors = [nbr.GetIdx() for nbr in mol.GetAtomWithIdx(ket_idx).GetNeighbors() if nbr.GetIdx() in ring_set]
            if dbond_a in ket_neighbors or dbond_b in ket_neighbors:
                return True, "Found a cyclohexenone: six‐membered non‐aromatic ring with one internal double bond and conjugated ketone"
    
    return False, "No cyclohexenone ring (six‑membered, non‑aromatic with one double bond and conjugated ketone) found"


# Example usage: testing with simple candidates (the provided examples are complex)
if __name__ == "__main__":
    test_smiles = [
        "O=C1CCCC=C1",       # cyclohex-2-enone, should be True
        "O=C1CC=CCC1",       # alternative connectivity, should be True
        "O=C1C=CCCC1",       # slight variation in connectivity, should be True; note that not all valid depictions are captured
    ]
    for smi in test_smiles:
        classified, reason = is_cyclohexenones(smi)
        print(f"SMILES: {smi}\nIs cyclohexenone? {classified}\nReason: {reason}\n")