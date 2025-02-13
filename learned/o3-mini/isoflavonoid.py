"""
Classifies: CHEBI:50753 isoflavonoid
"""
"""
Classifies: isoflavonoid
Definition: Any 1-benzopyran with an aryl substituent at position 3.
The algorithm first obtains the ring information and looks for two fused rings that meet the following:
  – One ring is a benzene ring (6 atoms, all carbons, aromatic).
  – The other is a “pyran” ring (6 atoms, exactly one oxygen and five carbons).
The two rings must share exactly two atoms (i.e. be fused).
Then the algorithm looks at the “free” (non‐fused) atoms of the pyran ring and, for each free atom that is a carbon,
checks if it has a neighbor (via a single bond) that is aromatic and is itself part of a benzene ring (i.e. an aryl substituent).
If so, it classifies the molecule as an isoflavonoid.
Note: This heuristic may not capture every isoflavonoid but improves on previous false positives/negatives.
"""

from rdkit import Chem

def is_isoflavonoid(smiles: str):
    """
    Determines if a molecule is an isoflavonoid based on its SMILES string.
    An isoflavonoid (for our purpose) is defined as a molecule containing a 1-benzopyran core --
    that is, a fused bicyclic system consisting of a benzene ring (6-membered aromatic,
    all carbons) fused (sharing exactly 2 atoms) with a pyran ring (6-membered containing 1 oxygen and 5 carbons)
        with an aryl substituent (an extra aromatic ring, typically a benzene ring) attached at a free (non-fused)
        position on the pyran (i.e. likely at position 3).
        
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as an isoflavonoid, False otherwise.
        str: Explanation of the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    ring_info = mol.GetRingInfo().AtomRings()
    if not ring_info:
        return False, "No rings found in molecule"
    
    # Helper functions for ring properties
    def is_benzene(ring):
        # A benzene ring: 6 atoms, all carbons that are aromatic.
        if len(ring) != 6:
            return False
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() != 6 or not atom.GetIsAromatic():
                return False
        return True

    def is_pyran(ring):
        # A pyran ring (for this purpose): 6 atoms with exactly 1 oxygen and 5 carbons.
        if len(ring) != 6:
            return False
        oxygen_count = 0
        carbon_count = 0
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() == 8:
                oxygen_count += 1
            elif atom.GetAtomicNum() == 6:
                carbon_count += 1
            else:
                return False
        return (oxygen_count == 1 and carbon_count == 5)
    
    # Look for two fused rings (sharing exactly 2 atoms) that meet our benzene and pyran criteria.
    for i in range(len(ring_info)):
        for j in range(i+1, len(ring_info)):
            ring1 = ring_info[i]
            ring2 = ring_info[j]
            shared = set(ring1).intersection(set(ring2))
            if len(shared) != 2:
                continue  # Skip if not fused by exactly 2 atoms.
            candidate = None
            if is_benzene(ring1) and is_pyran(ring2):
                benzene_ring = set(ring1)
                pyran_ring = set(ring2)
                candidate = True
            elif is_benzene(ring2) and is_pyran(ring1):
                benzene_ring = set(ring2)
                pyran_ring = set(ring1)
                candidate = True
            else:
                continue
            if candidate:
                # Identify the free (non-fused) atoms of the pyran ring.
                free_atoms = [idx for idx in pyran_ring if idx not in shared]
                # Look for at least one free pyran atom (expected at position 3) that has an aryl substituent
                # (a neighbor outside the pyran that is aromatic and belongs to a benzene ring).
                for idx in free_atoms:
                    atom = mol.GetAtomWithIdx(idx)
                    # Consider only carbon atoms in the pyran (the key substitution should be on a carbon).
                    if atom.GetAtomicNum() != 6:
                        continue
                    for nbr in atom.GetNeighbors():
                        if nbr.GetIdx() in pyran_ring:
                            continue  # Skip neighbors belonging to the pyran system.
                        # Check bond type if needed, here we only consider any connection.
                        if not nbr.GetIsAromatic():
                            continue
                        # Now verify that this neighbor is part of a benzene ring.
                        for ring in ring_info:
                            if nbr.GetIdx() in ring and len(ring) == 6 and all(mol.GetAtomWithIdx(a).GetAtomicNum()==6 and mol.GetAtomWithIdx(a).GetIsAromatic() for a in ring):
                                return True, ("Found benzopyran core (benzene fused with a 6-membered ring containing one oxygen) "
                                               "with an aryl substituent attached at a free pyran carbon (likely position 3).")
    return False, "No benzopyran core with an aryl substituent at position 3 was detected"

# Example usage (you can test with your compound SMILES):
if __name__ == "__main__":
    # Testing on daidzein (an isoflavonoid) and a negative control.
    test_smiles = [
        "Oc1ccc(cc1)-c1coc2cc(O)ccc2c1=O",  # daidzein (should be True)
        "CC(=O)O",  # Acetic acid (negative control)
    ]
    for smi in test_smiles:
        classification, reason = is_isoflavonoid(smi)
        print(f"SMILES: {smi}\nClassified as isoflavonoid: {classification}\nReason: {reason}\n")