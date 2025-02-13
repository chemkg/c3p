"""
Classifies: CHEBI:23763 pyrroline
"""
"""
Classifies: any organic heteromonocyclic compound based on a dihydropyrrole (pyrroline)
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_pyrroline(smiles: str):
    """
    Determines if a molecule is a pyrroline derivative (i.e. based on a dihydropyrrole ring)
    by checking for a five-membered ring that contains exactly one nitrogen and four carbons,
    and which has exactly one double bond.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if a dihydropyrrole ring is found, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."

    # Check that the molecule is organic (has at least one carbon atom)
    if not any(atom.GetAtomicNum() == 6 for atom in mol.GetAtoms()):
        return False, "Molecule lacks carbon atoms; not organic."

    # Get ring information from the molecule.
    ring_info = mol.GetRingInfo().AtomRings()
    
    # Iterate through each ring in the molecule.
    for ring in ring_info:
        if len(ring) == 5:  # Looking for five-membered rings
            # Count the number of nitrogen and carbon atoms in the ring.
            n_count = 0
            c_count = 0
            for idx in ring:
                atom = mol.GetAtomWithIdx(idx)
                if atom.GetAtomicNum() == 7:
                    n_count += 1
                elif atom.GetAtomicNum() == 6:
                    c_count += 1
            # Check that the ring has exactly one nitrogen and four carbons.
            if n_count == 1 and c_count == 4:
                # Now check the bonds within the ring for double bonds.
                double_bond_count = 0
                # Loop over each bond that connects consecutive atoms in the ring.
                for i in range(len(ring)):
                    a1 = ring[i]
                    a2 = ring[(i + 1) % len(ring)]  # wrap around for the last bond in the ring
                    bond = mol.GetBondBetweenAtoms(a1, a2)
                    if bond is None:
                        continue
                    if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                        double_bond_count += 1
                # We require exactly one double bond in the ring (dihydropyrrole motif).
                if double_bond_count == 1:
                    return True, "Found a five-membered dihydropyrrole ring (1 nitrogen, 4 carbons, 1 double bond)."
                    
    # If we reach here, no ring satisfied our criteria.
    return False, "No five-membered dihydropyrrole ring (1 nitrogen, 4 carbons, and 1 double bond) found."