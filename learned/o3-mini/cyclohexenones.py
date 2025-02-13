"""
Classifies: CHEBI:48953 cyclohexenones
"""
"""
Classifies: Cyclohexenones – Any six-membered alicyclic ketone having one double bond in the ring.
"""

from rdkit import Chem

def is_cyclohexenones(smiles: str):
    """
    Determines if a molecule is a cyclohexenone based on its SMILES string.
    A cyclohexenone is defined as a six-membered alicyclic (i.e. non-aromatic) ring, all carbon atoms,
    with one ketone functional group (a carbon (in the ring) double-bonded to O) and exactly one double bond
    (C=C) within the ring.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule meets criteria, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    ri = mol.GetRingInfo()
    ring_atom_indices = ri.AtomRings()
    
    # Iterate over rings to search for a matching six-membered ring.
    for ring in ring_atom_indices:
        if len(ring) != 6:
            # Only interested in rings of exactly six atoms.
            continue
        
        # Check that all atoms in the ring are carbons and none are aromatic.
        all_carbons = True
        non_aromatic = True
        for atom_idx in ring:
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetAtomicNum() != 6:
                all_carbons = False
                break
            if atom.GetIsAromatic():
                non_aromatic = False
                break
        if not all_carbons:
            continue
        if not non_aromatic:
            continue
        
        # Flag to check if the ketone (carbonyl) group is present on the ring:
        ketone_found = False
        
        # Count how many double bonds are present between ring atoms.
        ring_double_bond_count = 0
        
        # Create a set for quick look-up to check if a neighbor is in-ring.
        ring_set = set(ring)
        
        # Loop over the bonds that connect atoms inside the ring.
        # Also, check for the ketone functionality.
        for atom_idx in ring:
            atom = mol.GetAtomWithIdx(atom_idx)
            # Check neighbors for a carbonyl group (C=O) with the oxygen outside the ring.
            if atom.GetAtomicNum() == 6:
                for nbr in atom.GetNeighbors():
                    # If neighbor is not in the ring and is oxygen with a double bond then it qualifies.
                    if nbr.GetAtomicNum() == 8 and nbr.GetIdx() not in ring:
                        # Determine the bond between atom and nbr
                        bond = mol.GetBondBetweenAtoms(atom_idx, nbr.GetIdx())
                        if bond is not None and bond.GetBondTypeAsDouble() == 2:
                            ketone_found = True
            # Now iterate over bonds of this atom to count double bonds inside the ring.
            for bond in atom.GetBonds():
                # We only count the bond if both atoms are in this ring.
                begin = bond.GetBeginAtomIdx()
                end   = bond.GetEndAtomIdx()
                if begin in ring_set and end in ring_set:
                    # Count as double bond if bond order is 2.
                    if bond.GetBondTypeAsDouble() == 2:
                        ring_double_bond_count += 1
        # In a ring, each bond is encountered twice (once for each atom), so halve the count.
        ring_double_bond_count //= 2
        
        # Check if exactly one double bond in ring AND a ketone group is present.
        if ring_double_bond_count == 1 and ketone_found:
            return True, "Found a six-membered non-aromatic ring with one double bond and a ketone group attached"
    
    return False, "No six-membered alicyclic cyclohexenone ring found"
    
# Example usage:
if __name__ == "__main__":
    # Here are some test SMILES (the simplest being cyclohex-2-enone)
    smiles_examples = [
        "O=C1CCCC=C1",  # cyclohex-2-enone (should be True)
        "O=C1CC=CCC1",  # alternative depiction
        "O=C1C=CC2CCCCC2"  # not a simple cyclohexenone ring
    ]
    for smi in smiles_examples:
        result, reason = is_cyclohexenones(smi)
        print(f"SMILES: {smi}\nClassified as cyclohexenone: {result}, Reason: {reason}\n")