"""
Classifies: CHEBI:140331 4'-hydroxyflavanones
"""
"""
Classifies: 4'-hydroxyflavanones

Definition: Any hydroxyflavanone having a hydroxy substituent located at position 4' 
(i.e. on the para position of the B ring, relative to the point where the B ring attaches to the flavanone core).
The algorithm first locates a flavanone core (2-phenylchroman-4-one) and then searches for an 
external –OH group at the para position of the B ring.
"""

from rdkit import Chem

def is_4__hydroxyflavanones(smiles: str):
    """
    Determines if a molecule is a 4'-hydroxyflavanone based on its SMILES string.
    
    A 4'-hydroxyflavanone is defined as a molecule containing a flavanone core
    (2-phenylchroman-4-one) with an aromatic B ring having a free hydroxy (-OH) substituent 
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
    
    # 1. Look for the flavanone core: an approximation of the 2-phenylchroman-4-one scaffold.
    # This SMARTS expects a six-membered heterocycle with a ketone (C=O) fused to a benzene ring.
    core_pattern = Chem.MolFromSmarts("C1CC(=O)c2ccccc2O1")
    core_matches = mol.GetSubstructMatches(core_pattern)
    if not core_matches:
        return False, "Flavanone core (2-phenylchroman-4-one) not found"
    
    # Work with the atoms of the first core match.
    core_atoms = set(core_matches[0])
    
    # Get ring information from the molecule.
    ring_info = mol.GetRingInfo().AtomRings()
    
    # 2. Search for the bond connecting the core to the aromatic B ring.
    # We expect this bond to be between a core (sp³, non‐aromatic) atom and an aromatic atom (which is not in the core).
    found_para_hydroxy = False
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        idx1 = a1.GetIdx()
        idx2 = a2.GetIdx()
        
        # Check if one atom is in the core and nonaromatic while the other is aromatic and not in the core.
        if ((idx1 in core_atoms and (not a1.GetIsAromatic()) and a2.GetIsAromatic() and (idx2 not in core_atoms)) or
            (idx2 in core_atoms and (not a2.GetIsAromatic()) and a1.GetIsAromatic() and (idx1 not in core_atoms))):
            
            # Determine which is the core atom (expected sp³) and which is the aromatic attachment from the B ring.
            if idx1 in core_atoms:
                core_idx = idx1
                arom_idx = idx2
            else:
                core_idx = idx2
                arom_idx = idx1
            
            # 3. Identify the aromatic ring (the B ring) that contains the aromatic attachment.
            for ring in ring_info:
                if len(ring) != 6:
                    continue  # Only interested in six-membered rings (benzene)
                # Ensure the ring is fully aromatic.
                if not all(mol.GetAtomWithIdx(ai).GetIsAromatic() for ai in ring):
                    continue
                if arom_idx not in ring:
                    continue
                
                # In a benzene ring, the substituent at the para position is three bonds away in the ring.
                pos = ring.index(arom_idx)
                para_idx = ring[(pos + 3) % 6]
                para_atom = mol.GetAtomWithIdx(para_idx)
                
                # 4. Check if the para_atom carries a free hydroxy group.
                # Look for a neighboring oxygen atom that is attached externally (i.e. not in the ring)
                # and that qualifies as a -OH (i.e. the oxygen should have only one heavy-atom neighbor).
                for nb in para_atom.GetNeighbors():
                    if nb.GetIdx() in ring:
                        continue  # Ignore atoms already part of the ring
                    if nb.GetAtomicNum() == 8:
                        # Check if the oxygen looks like a free -OH.
                        # For a free hydroxy, the oxygen should have only one neighbor (the ring carbon)
                        # and at least one implicit hydrogen.
                        if nb.GetDegree() == 1 and nb.GetTotalNumHs() >= 1:
                            found_para_hydroxy = True
                            break
                if found_para_hydroxy:
                    break
            if found_para_hydroxy:
                break

    if not found_para_hydroxy:
        return False, "No free hydroxy group found at the para (4') position on the B ring"
    
    return True, "Contains flavanone core with a free -OH substituent at the 4' (para) position on the B ring"


# Optional testing code (can be removed or modified as needed).
if __name__ == '__main__':
    test_smiles = [
        "CC(C)=CCc1c(O)ccc2C(=O)CC(Oc12)c1ccc(O)cc1O",  # euchrenone-a7 (True)
        "O[C@@H]1[C@H](Oc2cc(O)cc(O)c2C1=O)c1ccc(O)c(O)c1",  # (+)-taxifolin (True)
        "Oc1ccc(cc1)[C@@H]1CC(=O)c2ccccc2O1",  # (2S)-4'-hydroxyflavanone (True)
        "O1[C@H](CC(=O)C=2C1=CC(O)=CC2O)C3=CC(O)=C(OC)C=C3"  # false positive candidate from previous attempt
    ]
    for s in test_smiles:
        result, reason = is_4__hydroxyflavanones(s)
        print(f"SMILES: {s}\nResult: {result}, Reason: {reason}\n")