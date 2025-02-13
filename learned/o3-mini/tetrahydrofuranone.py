"""
Classifies: CHEBI:47016 tetrahydrofuranone
"""
"""
Classifies: tetrahydrofuranone
Definition: Any oxolane (i.e. a 5-membered saturated ring containing one oxygen and four carbons)
with an oxo- substituent (a carbonyl group) attached on one of the ring carbons.
The oxo substituent may be attached directly (via a C=O) or through an acyl group.
We now ensure that the candidate oxolane ring is fully saturated (non‐aromatic, single bonds only)
to avoid false positives from aromatic or fused rings.
"""

from rdkit import Chem
from rdkit.Chem import rdchem

def is_tetrahydrofuranone(smiles: str):
    """
    Determines if a molecule is a tetrahydrofuranone based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        (bool, str): A tuple where the boolean is True if the molecule is classified as a tetrahydrofuranone,
                     and the string provides the reason for the classification.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get ring information.
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    
    # Iterate over candidate rings.
    for ring in rings:
        # Only consider 5-membered rings.
        if len(ring) != 5:
            continue

        # Count the number of oxygen and carbon atoms in the candidate ring.
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

        # Check every bond between adjacent atoms in the ring.
        ring_bonds_saturated = True
        n = len(ring)
        for i in range(n):
            a1 = ring[i]
            a2 = ring[(i+1) % n]  # wrap-around for cyclicity
            bond = mol.GetBondBetweenAtoms(a1, a2)
            if bond is None or bond.GetBondType() != rdchem.BondType.SINGLE:
                ring_bonds_saturated = False
                break
        if not ring_bonds_saturated:
            continue

        # Now, look for an oxo substituent on one of the ring carbons.
        # We check two possibilities:
        #  1. A ring carbon directly bears an exocyclic oxygen attached by a double bond (a carbonyl).
        #  2. A ring carbon is attached to a carbon (acyl group) that carries a double bonded oxygen.
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
                # Possibility #1: Direct exocyclic carbonyl group.
                if nbr.GetAtomicNum() == 8 and bond.GetBondType() == rdchem.BondType.DOUBLE:
                    return True, "Found a tetrahydrofuran ring with a direct carbonyl substituent on a ring carbon"
                # Possibility #2: Acyl substituent attached via a carbon.
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
    return False, "No tetrahydrofuran ring with an oxo substituent found"

# Example usage (you may add additional tests here):
if __name__ == "__main__":
    # Test example: N-isovaleryl-L-homoserine lactone
    test_smiles = "CC(C)CC(=O)N[C@H]1CCOC1=O"
    result, reason = is_tetrahydrofuranone(test_smiles)
    print(f"SMILES: {test_smiles}")
    print("Classification:", result)
    print("Reason:", reason)