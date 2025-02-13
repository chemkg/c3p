"""
Classifies: CHEBI:47016 tetrahydrofuranone
"""
"""
Classifies: tetrahydrofuranone
Definition: Any oxolane (i.e. a 5‐membered saturated ring with one oxygen and four carbons)
that carries an oxo- substituent (a carbonyl, either as a direct C=O or through an acyl group)
on one of the ring carbons. Additionally, the oxolane ring must be isolated (i.e. not fused with other rings)
to avoid false positives.
"""

from rdkit import Chem
from rdkit.Chem import rdchem

def is_tetrahydrofuranone(smiles: str):
    """
    Determines if a molecule qualifies as a tetrahydrofuranone based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        (bool, str): Tuple where the boolean indicates if the molecule is classified as a tetrahydrofuranone,
                     and the string provides the reason for the classification.
    """
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Retrieve ring information.
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()

    # Iterate over rings to find a candidate oxolane ring.
    for ring in rings:
        # Only consider 5-membered rings.
        if len(ring) != 5:
            continue

        # Count atoms in the candidate ring: need exactly four carbons and one oxygen.
        oxygen_count = 0
        carbon_count = 0
        for idx in ring:
            atomic_num = mol.GetAtomWithIdx(idx).GetAtomicNum()
            if atomic_num == 8:
                oxygen_count += 1
            elif atomic_num == 6:
                carbon_count += 1
        if oxygen_count != 1 or carbon_count != 4:
            continue

        # Ensure the ring is saturated:
        # a) none of its atoms is aromatic,
        # b) all bonds between atoms in the ring are single bonds.
        saturated = True
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetIsAromatic():
                saturated = False
                break
        if not saturated:
            continue

        ring_bonds_saturated = True
        n = len(ring)
        for i in range(n):
            a1 = ring[i]
            a2 = ring[(i + 1) % n]  # wrap around for cyclicity
            bond = mol.GetBondBetweenAtoms(a1, a2)
            if bond is None or bond.GetBondType() != rdchem.BondType.SINGLE:
                ring_bonds_saturated = False
                break
        if not ring_bonds_saturated:
            continue

        # NEW CHECK: Ensure that the candidate oxolane ring is isolated.
        # Every atom in the ring should belong to exactly one ring.
        isolated = True
        for idx in ring:
            if ring_info.NumAtomRings(idx) != 1:
                isolated = False
                break
        if not isolated:
            continue

        # Now check for an oxo substituent on one of the ring carbons.
        # Possibility 1: The ring carbon directly bears an exocyclic oxygen by a double bond.
        # Possibility 2: The ring carbon is attached to a carbon (acyl group) that carries a C=O.
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            # Only consider carbon atoms in the ring.
            if atom.GetAtomicNum() != 6:
                continue

            # Inspect neighbors that are not part of the ring.
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring:
                    continue
                bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
                if bond is None:
                    continue
                # Possibility 1: Direct carbonyl group (C=O directly attached to ring carbon).
                if nbr.GetAtomicNum() == 8 and bond.GetBondType() == rdchem.BondType.DOUBLE:
                    return True, "Found a tetrahydrofuran ring with a direct carbonyl substituent on a ring carbon"
                # Possibility 2: Acyl substituent. Check that the neighbor is a carbon
                # that carries a double bonded oxygen.
                if nbr.GetAtomicNum() == 6:
                    for sub_nbr in nbr.GetNeighbors():
                        if sub_nbr.GetIdx() == atom.GetIdx():
                            continue
                        sub_bond = mol.GetBondBetweenAtoms(nbr.GetIdx(), sub_nbr.GetIdx())
                        if sub_bond is None:
                            continue
                        if sub_nbr.GetAtomicNum() == 8 and sub_bond.GetBondType() == rdchem.BondType.DOUBLE:
                            return True, "Found a tetrahydrofuran ring with an acyl substituent on a ring carbon"
    # If no candidate ring passed all checks, classification fails.
    return False, "No isolated tetrahydrofuran ring with an oxo substituent found"

# Example usage (feel free to add additional tests):
if __name__ == "__main__":
    # Example: N-isovaleryl-L-homoserine lactone
    test_smiles = "CC(C)CC(=O)N[C@H]1CCOC1=O"
    result, reason = is_tetrahydrofuranone(test_smiles)
    print(f"SMILES: {test_smiles}")
    print("Classification:", result)
    print("Reason:", reason)