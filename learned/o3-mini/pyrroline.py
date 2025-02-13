"""
Classifies: CHEBI:23763 pyrroline
"""
"""
Classifies: any organic heteromonocyclic compound with a structure based on a dihydropyrrole (pyrroline)

Definition: A pyrroline derivative (dihydropyrrole) is defined as a five-membered ring containing exactly one nitrogen 
and four carbons with exactly one double bond among the ring bonds.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_pyrroline(smiles: str):
    """
    Determines if a molecule is a pyrroline derivative (i.e. based on a dihydropyrrole ring)
    by checking for a five-membered ring that contains exactly one nitrogen and four carbons,
    with exactly one double bond in that ring.

    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if a valid dihydropyrrole ring is found, False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
        
    # Ensure the molecule is organic (has at least one carbon atom).
    if not any(atom.GetAtomicNum() == 6 for atom in mol.GetAtoms()):
        return False, "Molecule lacks carbon atoms; not organic."
        
    # Force Kekulization to get explicit bond orders.
    try:
        Chem.Kekulize(mol, clearAromaticFlags=True)
    except Exception as e:
        return False, "Kekulization failed; ambiguous aromaticity."
    
    # Retrieve ring information: both atom rings and corresponding bond rings.
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()  # tuple of tuples with atom indices
    bond_rings = ring_info.BondRings()  # tuple of tuples with bond indices corresponding to the atom rings

    # Iterate over each ring candidate.
    for i, ring in enumerate(atom_rings):
        # Consider only five-membered rings.
        if len(ring) != 5:
            continue

        # Count the number of nitrogen and carbon atoms in the ring.
        n_count = 0
        c_count = 0
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            atomic_num = atom.GetAtomicNum()
            if atomic_num == 7:
                n_count += 1
            elif atomic_num == 6:
                c_count += 1
        # Skip if not exactly 1 nitrogen and 4 carbons.
        if n_count != 1 or c_count != 4:
            continue

        # Get the corresponding bond indices for the ring.
        # bond_rings[i] returns the bonds that complete this specific ring.
        bonds_in_ring = bond_rings[i]

        # Count the number of double bonds among those bonds.
        double_bond_count = 0
        for bond_idx in bonds_in_ring:
            bond = mol.GetBondWithIdx(bond_idx)
            if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                double_bond_count += 1
        
        # Check if the ring contains exactly one double bond.
        if double_bond_count == 1:
            return True, ("Found a valid five-membered dihydropyrrole ring with "
                          "1 nitrogen, 4 carbons, and exactly 1 double bond.")
    
    return False, ("No five-membered dihydropyrrole ring found that has exactly "
                   "1 nitrogen, 4 carbons, and 1 double bond.")

# Example usage for internal testing.
if __name__ == "__main__":
    test_smiles = [
        "O=C/1NCC(\\C1=C(\\O)/C=C/C(=C/C#C/C=C/C)/C)=O",  # Ravynic acid (expected: True)
        "S=C1NCCC1",                                    # Pyrrolidine-2-thione (expected: False: no double bond)
        "C1(N=CCC1)(C)C",                               # 5,5-dimethyl-1-pyrroline (expected: True)
        "OC1=C2[N+](=CC=C1)[C@@H](/C=C(/C=C/C3CC3)\\C)[C@@H]([C@]2(O)C)O",  # Cyclizidine F (expected: True)
        "O=C1N(C(O)(CC)C(=C1C(=O)C(CCCCCCC)C)O)C",      # Penicillenol D (expected: True)
        "O=C1N[C@](C(=O)N)(CC(C)C)[C@](C1)(O)C",         # Monascustin (expected: False)
    ]
    
    for smi in test_smiles:
        result, reason = is_pyrroline(smi)
        print(f"SMILES: {smi}\nResult: {result}\nReason: {reason}\n{'-'*40}")