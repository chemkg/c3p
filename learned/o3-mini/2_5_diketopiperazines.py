"""
Classifies: CHEBI:65061 2,5-diketopiperazines
"""
"""
Classifies: CHEBI:2,5-diketopiperazines
Definition: Any piperazinone that has a piperazine-2,5-dione skeleton.
This implementation searches for a six-membered ring (from the molecule's symmetric SSSR)
in which all atoms are either carbon or nitrogen, with exactly two nitrogen atoms and four
carbon atoms. It then checks that exactly two of the ring carbons are carbonyl (C=O)
(with a double bond to oxygen) and that in the ring ordering these two carbonyl carbons
are separated by three bonds (i.e. they are in the 2 and 5 positions).
"""

from rdkit import Chem

def is_2_5_diketopiperazines(smiles: str):
    """
    Determines if a molecule contains a piperazine-2,5-dione core (a 2,5-diketopiperazine)
    by examining its ring systems.
    
    For each six-membered ring:
      - Checks that every atom is either a carbon (atomic number 6) or nitrogen (7).
      - Confirms that exactly 2 nitrogens and 4 carbons are present.
      - Identifies which ring carbons have a double-bonded oxygen (carbonyl).
      - Accepts a ring when exactly 2 carbonyl-bearing carbons are found, and their
        positions in the ring (via the cycle ordering) are separated by 3 (i.e. have two atoms
        between them on one side of the ring, corresponding to positions 2 and 5 in a 1-indexed cycle).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule contains a complete piperazine-2,5-dione ring.
        str: Reason for the classification decision.
    """
    
    # Convert SMILES to RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # Get all rings as determined by the symmetric SSSR algorithm.
    rings = Chem.GetSymmSSSR(mol)
    
    # Loop over each ring.
    for ring in rings:
        if len(ring) != 6:
            continue  # only consider 6-membered rings
        
        # The ring is a tuple of atom indices in a cyclic order.
        ring_atoms = list(ring)
        
        # Check that every atom in the ring is C (6) or N (7)
        atom_nums = [mol.GetAtomWithIdx(idx).GetAtomicNum() for idx in ring_atoms]
        if any(num not in (6, 7) for num in atom_nums):
            continue

        # Verify that exactly two nitrogens and four carbons are present.
        if atom_nums.count(7) != 2 or atom_nums.count(6) != 4:
            continue
        
        # Identify positions in the ring (by order) of carbons that are carbonyl (C=O).
        carbonyl_positions = []
        for pos, idx in enumerate(ring_atoms):
            atom = mol.GetAtomWithIdx(idx)
            # Only consider carbon atoms.
            if atom.GetAtomicNum() != 6:
                continue
            # Search bonds for a double bond to oxygen.
            for bond in atom.GetBonds():
                # Check bond type is double (a carbonyl bond)
                if bond.GetBondType() == Chem.BondType.DOUBLE:
                    neighbor = bond.GetOtherAtom(atom)
                    if neighbor.GetAtomicNum() == 8:
                        carbonyl_positions.append(pos)
                        break
        # We require exactly two carbonyls on ring carbons.
        if len(carbonyl_positions) != 2:
            continue
        
        # Check that the two carbonyl positions are spaced three bonds apart in the ring.
        pos1, pos2 = carbonyl_positions
        diff = abs(pos1 - pos2)
        cyclic_distance = min(diff, 6 - diff)
        if cyclic_distance == 3:
            return True, "Molecule contains a complete piperazine-2,5-dione skeleton."
    
    return False, "Molecule does not contain a complete piperazine-2,5-dione ring."


# Example usage (you can add more tests as needed)
if __name__ == "__main__":
    test_molecules = [
        ("Brocazine F", "S1S[C@]23N([C@@H]4[C@@H](O)C=C[C@@H]([C@H]4C2)O)C([C@]15N([C@@H]6[C@@H](O)C=CC([C@H]6C5)=O)C3=O)=O"),
        ("mycocyclosin", "Oc1ccc2C[C@@H]3NC(=O)[C@H](Cc4ccc(O)c(c4)-c1c2)NC3=O"),
        ("tardioxopiperazine A", "O=C1N[C@H](C(=O)N[C@H]1CC=2C=3C(NC2C(C=C)(C)C)=CC=C(C3)CC=C(C)C)C"),
        ("piperazine-2,5-dione", "O=C1CNC(=O)CN1")
    ]
    
    for name, smi in test_molecules:
        result, reason = is_2_5_diketopiperazines(smi)
        print(f"{name}:\n  Result: {result}\n  Reason: {reason}\n")