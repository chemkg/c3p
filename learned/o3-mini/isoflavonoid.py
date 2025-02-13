"""
Classifies: CHEBI:50753 isoflavonoid
"""
"""
Classifies: isoflavonoid
Definition: Any 1-benzopyran with an aryl substituent at position 3.
A 1-benzopyran is defined here as a bicyclic system comprising a pyran ring 
(6-membered: 5 carbons and 1 oxygen) fused to a benzene ring (6 aromatic carbons)
via exactly 2 atoms. In typical isoflavonoid numbering the oxygen is at position 1,
the fused atom that is bonded to that oxygen is position 2 and the other fused atom 
(which is not bonded directly to the oxygen) is position 3. An isoflavonoid further 
requires that the position 3 atom carries one—and only one—extra substituent that is 
an aryl group (that is, a benzene ring, as defined by 6 aromatic carbons).
Note: This heuristic does not cover every edge case.
"""

from rdkit import Chem

def is_isoflavonoid(smiles: str):
    """
    Determines if a molecule is an isoflavonoid according to our definition.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if classified as isoflavonoid, False otherwise.
        str: Explanation of the decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    ring_info = mol.GetRingInfo().AtomRings()
    if not ring_info:
        return False, "No rings found in molecule"
    
    # Helper: check if a ring is benzene-like: exactly 6 atoms, all aromatic carbons.
    def is_benzene(ring):
        if len(ring) != 6:
            return False
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() != 6 or not atom.GetIsAromatic():
                return False
        return True

    # Helper: check if a ring is pyran-like:
    # exactly 6 atoms, five carbons and one oxygen.
    def is_pyran(ring):
        if len(ring) != 6:
            return False
        o_count = 0
        c_count = 0
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() == 8:
                o_count += 1
            elif atom.GetAtomicNum() == 6:
                c_count += 1
            else:
                return False
        return (o_count == 1 and c_count == 5)
    
    # Helper: check if an atom (by index) belongs to an aromatic benzene ring.
    def atom_in_benzene(idx):
        for ring in ring_info:
            if len(ring) == 6 and idx in ring and is_benzene(ring):
                return True
        return False

    # We now iterate over each pair of rings in the molecule.
    # We look for a fused pair (sharing exactly 2 atoms) where one ring is benzene
    # and the other is pyran.
    for i in range(len(ring_info)):
        for j in range(i+1, len(ring_info)):
            ring1 = set(ring_info[i])
            ring2 = set(ring_info[j])
            shared = ring1.intersection(ring2)
            if len(shared) != 2:
                continue  # not fused by exactly two atoms
            # Identify which ring is benzene and which is pyran.
            benzene_ring = None
            pyran_ring = None
            if is_benzene(ring_info[i]) and is_pyran(ring_info[j]):
                benzene_ring = set(ring_info[i])
                pyran_ring = set(ring_info[j])
            elif is_benzene(ring_info[j]) and is_pyran(ring_info[i]):
                benzene_ring = set(ring_info[j])
                pyran_ring = set(ring_info[i])
            else:
                continue
            
            # For our fused system, we expect the pyran ring to contain exactly one oxygen.
            o_atoms = [idx for idx in pyran_ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8]
            if len(o_atoms) != 1:
                continue
            o_idx = o_atoms[0]
            
            # Our two fused atoms come from the intersection.
            if len(shared) != 2:
                continue
            shared_list = list(shared)
            # Among the two fused atoms, the one directly bonded to the oxygen in the pyran ring
            # is considered to be position 2; the other will be candidate for position 3.
            pos2 = None
            pos3 = None
            for idx in shared_list:
                atom = mol.GetAtomWithIdx(idx)
                # Check if oxygen (o_idx) is a neighbor in the pyran ring.
                nbrs = [nbr.GetIdx() for nbr in atom.GetNeighbors()]
                if o_idx in nbrs:
                    pos2 = idx
                else:
                    pos3 = idx
            if pos3 is None:
                # Could not differentiate the two fused atoms.
                continue

            # Now check that the candidate position 3 atom carries one substituent 
            # (i.e. a neighbor that lies outside the pyran core) that is an aryl group.
            candidate_atom = mol.GetAtomWithIdx(pos3)
            substituent_count = 0
            for nbr in candidate_atom.GetNeighbors():
                nbr_idx = nbr.GetIdx()
                # Skip neighbors that are part of the pyran ring (the core).
                if nbr_idx in pyran_ring:
                    continue
                # We require that the bond is a single bond.
                bond = mol.GetBondBetweenAtoms(candidate_atom.GetIdx(), nbr_idx)
                if bond is None or bond.GetBondTypeAsDouble() != 1.0:
                    continue
                # Now, check that the neighbor is aromatic and belongs to a benzene ring.
                if not nbr.GetIsAromatic():
                    continue
                if atom_in_benzene(nbr_idx):
                    substituent_count += 1
            
            if substituent_count == 1:
                return True, ("Found benzopyran core (pyran fused to a benzene ring) with an extra aryl substituent "
                              "on the fused atom that is not bonded to the pyran oxygen (position 3).")
            elif substituent_count > 1:
                return False, (f"Multiple ({substituent_count}) aryl substituents found on the candidate position 3; "
                               "expected exactly one.")
    # If no valid fused benzopyran core with proper substitution was found:
    return False, "No valid benzopyran core with a 3-aryl substituent was detected"


# Example usage (for testing the function):
if __name__ == "__main__":
    test_smiles = [
        # These are examples of known or putative isoflavonoids:
        "Oc1ccc(cc1)-c1coc2cc(O)ccc2c1=O",  # daidzein (should be True)
        "O1C2=C(C=C(C)C)C(O)=CC(O)=C2C(=O)C(C3=CC=C(OC)C=C3)=C1",  # Gancaonin M (True)
        "O1C2=C(C/C=C(/CO)\\C)C(O)=CC(O)=C2C(=O)C(C3=C(O)C=C(O)C=C3)=C1",  # 2,3-Dehydrokievitol (True)
        # A negative control:
        "CC(=O)O",  # acetic acid (False)
    ]
    for smi in test_smiles:
        classification, reason = is_isoflavonoid(smi)
        print(f"SMILES: {smi}\nClassified as isoflavonoid: {classification}\nReason: {reason}\n")