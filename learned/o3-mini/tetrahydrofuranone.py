"""
Classifies: CHEBI:47016 tetrahydrofuranone
"""
"""
Classifies: tetrahydrofuranone
Definition: Any oxolane (five-membered ring comprising one oxygen and four carbons)
with an oxo- substituent (i.e. a carbonyl group) attached on one of the ring carbons.
We further require that the tetrahydrofuran ring be “isolated” (not fused to another 5-membered ring),
in an effort to avoid false positives from polycyclic or fused structures.
"""

from rdkit import Chem

def is_tetrahydrofuranone(smiles: str):
    """
    Determines if a molecule is a tetrahydrofuranone based on its SMILES string.
    
    A tetrahydrofuranone is defined here as any molecule that contains an oxolane
    ring (a 5-membered ring with exactly one oxygen and four carbons, not fused to another such ring)
    and such that at least one of the ring’s carbon atoms bears an oxo (-C=O) substituent.
    
    The oxo substituent is detected either directly (an exocyclic oxygen double-bonded to the ring carbon)
    or through an acyl group attached to the ring carbon.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a tetrahydrofuranone, False otherwise.
        str: Reason for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    ri = mol.GetRingInfo()
    rings = ri.AtomRings()
    
    # Iterate over candidate rings
    for ring in rings:
        if len(ring) != 5:
            continue  # Only consider 5-membered rings
        
        # Count atoms in the ring: we want exactly one oxygen and four carbons.
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
        
        # Check that the candidate 5-membered ring is isolated (not fused with another 5-membered ring)
        # Here, if any atom belongs to more than one 5-membered ring then we skip this candidate.
        is_fused = False
        for idx in ring:
            count_in_5 = sum(1 for r in rings if len(r)==5 and idx in r)
            if count_in_5 > 1:
                is_fused = True
                break
        if is_fused:
            continue
        
        # Now search for an oxo substituent on one of the ring carbons.
        # We check two possibilities:
        # 1. The ring carbon directly has an exocyclic oxygen attached via a double bond.
        # 2. The ring carbon is attached to an exocyclic carbon that itself carries a carbonyl group.
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            # Only check carbon atoms
            if atom.GetAtomicNum() != 6:
                continue
            for nbr in atom.GetNeighbors():
                # Consider only substituents not within the ring.
                if nbr.GetIdx() in ring:
                    continue
                bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
                if bond is None:
                    continue
                # Possibility #1: Direct exocyclic oxygen (C=O) attached
                if nbr.GetAtomicNum() == 8 and bond.GetBondTypeAsDouble() >= 2.0:
                    return True, "Found a tetrahydrofuran ring with a carbonyl substituent on a ring carbon"
                # Possibility #2: Exocyclic acyl substituent. The neighbor is a carbon that bears a carbonyl group.
                if nbr.GetAtomicNum() == 6:
                    for sub_nbr in nbr.GetNeighbors():
                        # Avoid backtracking to the ring carbon we came from.
                        if sub_nbr.GetIdx() == atom.GetIdx():
                            continue
                        sub_bond = mol.GetBondBetweenAtoms(nbr.GetIdx(), sub_nbr.GetIdx())
                        if sub_bond and sub_bond.GetBondTypeAsDouble() >= 2.0 and sub_nbr.GetAtomicNum() == 8:
                            return True, "Found a tetrahydrofuran ring with an acyl substituent on a ring carbon"
    return False, "No tetrahydrofuran ring with an oxo substituent found"

# Example usage (you may add additional tests here):
if __name__ == "__main__":
    # Example: N-isovaleryl-L-homoserine lactone         
    test_smiles = "CC(C)CC(=O)N[C@H]1CCOC1=O"
    result, reason = is_tetrahydrofuranone(test_smiles)
    print(f"SMILES: {test_smiles}")
    print("Classification:", result)
    print("Reason:", reason)