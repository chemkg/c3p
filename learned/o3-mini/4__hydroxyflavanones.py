"""
Classifies: CHEBI:140331 4'-hydroxyflavanones
"""
"""
Classifies: 4'-hydroxyflavanones
Definition: Any hydroxyflavanone having a hydroxy substituent located at position 4' (i.e. on the para position of the B ring).
This implementation first checks for the flavanone core (2-phenylchroman-4-one) using a SMARTS pattern.
Then, it searches for an aromatic ring (of 6 atoms) attached via a non‐aromatic (sp3) atom (the B ring)
and checks that the ring contains an –OH substituent in the para (3 bonds away) position relative to the point of attachment.
"""

from rdkit import Chem

def is_4__hydroxyflavanones(smiles: str):
    """
    Determines if a molecule is a 4'-hydroxyflavanone based on its SMILES string.
    
    A 4'-hydroxyflavanone is defined as a molecule containing a flavanone core (2-phenylchroman-4-one)
    with an aromatic B ring having an -OH group at the para position relative to its connection.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule qualifies as a 4'-hydroxyflavanone, False otherwise.
        str: Reason for the classification decision.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 1. Check for the flavanone core.
    # The core pattern "C1CC(=O)c2ccccc2O1" is an approximation of the 2-phenylchroman-4-one skeleton.
    core_pattern = Chem.MolFromSmarts("C1CC(=O)c2ccccc2O1")
    core_matches = mol.GetSubstructMatches(core_pattern)
    if not core_matches:
        return False, "Flavanone core (2-phenylchroman-4-one) not found"
    
    # 2. Look for the B ring connection and verify its para substituent.
    # We expect the B ring (a benzene ring) to be attached to the core via a non-aromatic (likely sp3) carbon.
    # We then check that within that benzene ring, the atom opposite (para) the attachment has an -OH.
    ring_info = mol.GetRingInfo().AtomRings()
    # Flag to record if we found the required substituent.
    found_para_hydroxy = False
    
    # For each bond, see if it connects a candidate core atom to an aromatic atom not in the core.
    # For simplicity, we use one of our core matches (they should share the same connectivity).
    core_atoms = set(core_matches[0])
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        idx1 = a1.GetIdx()
        idx2 = a2.GetIdx()
        # Identify candidate: one end is in the core and is non-aromatic (likely sp3) and the other is aromatic.
        if (idx1 in core_atoms and not a1.GetIsAromatic() and a2.GetIsAromatic() and idx2 not in core_atoms) or \
           (idx2 in core_atoms and not a2.GetIsAromatic() and a1.GetIsAromatic() and idx1 not in core_atoms):
            # Determine the aromatic atom (on the B ring) that is connected.
            if idx1 in core_atoms:
                core_idx = idx1
                arom_idx = idx2
            else:
                core_idx = idx2
                arom_idx = idx1
            
            # Now, look at rings that include arom_idx.
            for ring in ring_info:
                # Only consider six-membered rings.
                if len(ring) != 6:
                    continue
                # Confirm that every atom in the ring is aromatic.
                if not all(mol.GetAtomWithIdx(ai).GetIsAromatic() for ai in ring):
                    continue
                # Check if our aromatic attachment is in this ring.
                if arom_idx not in ring:
                    continue
                # Find the index of arom_idx in the ring list.
                pos = ring.index(arom_idx)
                # In a benzene ring, the atom para (3 bonds away) is at index (pos+3) mod 6.
                para_idx = ring[(pos + 3) % 6]
                para_atom = mol.GetAtomWithIdx(para_idx)
                # Check if para_atom has an -OH substituent.
                # We do this by examining its neighbors that are not part of the ring.
                for nb in para_atom.GetNeighbors():
                    # Skip if neighbor is also a ring atom (part of the aromatic system).
                    if nb.GetIdx() in ring:
                        continue
                    if nb.GetAtomicNum() == 8:  # oxygen
                        # Optionally, one could check bond type and degree.
                        # Here we assume that an -O attached (likely with an implicit H) acts as a hydroxy.
                        found_para_hydroxy = True
                        break
                if found_para_hydroxy:
                    break
            if found_para_hydroxy:
                break

    if not found_para_hydroxy:
        return False, "No hydroxyl group found at the para position on the B ring"
    
    return True, "Contains flavanone core with an -OH substituent at the 4' (para) position on the B ring"

# The function can be tested using the examples provided.
if __name__ == '__main__':
    test_smiles = [
        "CC(C)=CCc1c(O)ccc2C(=O)CC(Oc12)c1ccc(O)cc1O",  # euchrenone-a7 (should be True)
        "O[C@@H]1[C@H](Oc2cc(O)cc(O)c2C1=O)c1ccc(O)c(O)c1",  # (+)-taxifolin (should be True)
        "O1[C@H](CC(=O)C=2C1=CC(O)=CC2O)C3=CC(O)=C(OC)C=C3"  # A false positive candidate
    ]
    for s in test_smiles:
        result, reason = is_4__hydroxyflavanones(s)
        print(f"SMILES: {s}\nResult: {result}, Reason: {reason}\n")