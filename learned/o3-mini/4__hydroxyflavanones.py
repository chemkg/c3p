"""
Classifies: CHEBI:140331 4'-hydroxyflavanones
"""
"""
Classifies: 4'-hydroxyflavanones
Definition: Any hydroxyflavanone having a hydroxy substituent located at position 4'
(i.e. on the para position of the B ring).

This improved implementation attempts to minimize false positives by:
  1. Using a SMARTS that captures the flavanone core (2-phenylchroman-4-one).
  2. Selecting only bonds where a non‐aromatic (aliphatic, expected sp3)
     core carbon connects to an aromatic carbon that is not in the core.
  3. Identifying the specific six-membered aromatic ring of the B ring and checking
     that the atom para (three bonds away along the ring) to the point of connection
     carries a hydroxy substituent (–OH) attached outside the ring.
"""

from rdkit import Chem

def is_4__hydroxyflavanones(smiles: str):
    """
    Determines if a molecule is a 4'-hydroxyflavanone based on its SMILES string.
    
    A 4'-hydroxyflavanone is defined as a molecule containing a flavanone core
    (2-phenylchroman-4-one) with an aromatic B ring having a hydroxy (-OH) substituent 
    at the para (4') position relative to the point of attachment.
    
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
    # This SMARTS pattern is an approximation of the 2-phenylchroman-4-one scaffold.
    # It expects a six-membered oxygen containing ring with a ketone function (C=O),
    # fused to a benzene ring.
    core_pattern = Chem.MolFromSmarts("C1CC(=O)c2ccccc2O1")
    core_matches = mol.GetSubstructMatches(core_pattern)
    if not core_matches:
        return False, "Flavanone core (2-phenylchroman-4-one) not found"
    
    # For further steps, use the atoms from the first core match.
    core_atoms = set(core_matches[0])
    
    # 2. Look for the bond connecting the flavanone core to the B ring.
    # The expected connection is from an aliphatic (non‐aromatic) carbon of the core
    # (position 2, which is sp3) to an aromatic carbon (the attachment of the B ring).
    found_para_hydroxy = False
    ring_info = mol.GetRingInfo().AtomRings()
    
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        idx1 = a1.GetIdx()
        idx2 = a2.GetIdx()
        
        # Check for candidate bond:
        # One atom in the core must be non-aromatic (aliphatic) and the other not in core must be aromatic.
        if (idx1 in core_atoms and (not a1.GetIsAromatic()) and a2.GetIsAromatic() and (idx2 not in core_atoms)) or \
           (idx2 in core_atoms and (not a2.GetIsAromatic()) and a1.GetIsAromatic() and (idx1 not in core_atoms)):
            
            # Identify which atom is in the core (expected sp3) and which is the aromatic B-ring atom.
            if idx1 in core_atoms:
                core_idx = idx1
                arom_idx = idx2
            else:
                core_idx = idx2
                arom_idx = idx1

            # Now, find the aromatic ring that contains the aromatic attachment.
            # We expect the B ring to be a six-membered fully aromatic ring.
            for ring in ring_info:
                if len(ring) != 6:
                    continue  # Not a benzene ring
                # Check that every atom in the ring is aromatic.
                if not all(mol.GetAtomWithIdx(ai).GetIsAromatic() for ai in ring):
                    continue
                if arom_idx not in ring:
                    continue
                
                # Determine the position of the attachment atom within the ring.
                pos = ring.index(arom_idx)
                # In a benzene ring, the para substituent is located 3 bonds away (opposite in the ring)
                para_idx = ring[(pos + 3) % 6]
                para_atom = mol.GetAtomWithIdx(para_idx)
                
                # Check if the para_atom has a hydroxy substituent: that is,
                # one of its neighbors (not in the ring) is an oxygen.
                for nb in para_atom.GetNeighbors():
                    if nb.GetIdx() in ring:
                        continue  # skip ring atoms
                    if nb.GetAtomicNum() == 8:
                        # Found an oxygen attached externally
                        found_para_hydroxy = True
                        break
                if found_para_hydroxy:
                    break
            if found_para_hydroxy:
                break

    if not found_para_hydroxy:
        return False, "No hydroxy group found at the para (4') position on the B ring"
    
    return True, "Contains flavanone core with an -OH substituent at the 4' (para) position on the B ring"

# Example testing (these examples are provided in the task description).
if __name__ == '__main__':
    test_smiles = [
        "CC(C)=CCc1c(O)ccc2C(=O)CC(Oc12)c1ccc(O)cc1O",  # euchrenone-a7 (True)
        "O[C@@H]1[C@H](Oc2cc(O)cc(O)c2C1=O)c1ccc(O)c(O)c1",  # (+)-taxifolin (True)
        "Oc1ccc(cc1)[C@@H]1CC(=O)c2ccccc2O1",  # (2S)-4'-hydroxyflavanone (True)
        "O1[C@H](CC(=O)C=2C1=CC(O)=CC2O)C3=CC(O)=C(OC)C=C3"  # false positive candidate from before
    ]
    for s in test_smiles:
        result, reason = is_4__hydroxyflavanones(s)
        print(f"SMILES: {s}\nResult: {result}, Reason: {reason}\n")