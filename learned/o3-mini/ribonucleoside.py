"""
Classifies: CHEBI:18254 ribonucleoside
"""
"""
Classifies: Any nucleoside where the sugar component is D-ribose.
A ribonucleoside is defined as a nucleoside where the sugar is D-ribose – a five-membered (furanose) ring
composed of 4 carbons and 1 oxygen, with hydroxyl substituents (including an exocyclic CH2OH group) and a glycosidic bond
(from one of the ring carbons to a nucleobase nitrogen).
This algorithm extracts five-membered rings with exactly four carbons and one oxygen,
checks that at least one ring carbon is attached to a nitrogen (the nucleobase),
rejects rings that bear a phosphate (P), and finally verifies that the ring bears at least two alcohol groups
(and at least one CH2OH substituent) matching a ribose-like pattern.
"""
from rdkit import Chem

def is_ribonucleoside(smiles: str):
    """
    Determines if the molecule (given by SMILES string) is a ribonucleoside.
    The sugar portion must be D-ribose: a five-membered ring (4 C and 1 O) with appropriate hydroxyl substituents,
    one ring carbon attached to a CH2OH group, and one ring carbon linked via a C–N bond to a nucleobase.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if it is classified as a ribonucleoside, False otherwise.
        str: Reason for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens so we can count OH groups accurately.
    mol = Chem.AddHs(mol)
    
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    candidate_found = False
    
    # iterate over rings
    for ring in atom_rings:
        # We focus on rings of exactly 5 atoms.
        if len(ring) != 5:
            continue
        
        # Count atoms by type inside the ring.
        oxygen_count = 0
        carbon_count = 0
        for idx in ring:
            atomic_num = mol.GetAtomWithIdx(idx).GetAtomicNum()
            if atomic_num == 8:
                oxygen_count += 1
            elif atomic_num == 6:
                carbon_count += 1
        if oxygen_count != 1 or carbon_count != 4:
            continue  # not a ribofuranose-type ring
        
        # Check that none of the ring atoms is attached to phosphorus (exclude phosphate groups).
        phosphate_attached = False
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring:
                    continue
                if nbr.GetAtomicNum() == 15:
                    phosphate_attached = True
                    break
            if phosphate_attached:
                break
        if phosphate_attached:
            continue
        
        # Check for nucleobase: at least one ring carbon should be attached to a nitrogen outside the ring.
        nucleobase_attached = False
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() != 6:
                continue  # consider only carbons in the ring
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring:
                    continue
                if nbr.GetAtomicNum() == 7:
                    nucleobase_attached = True
                    break
            if nucleobase_attached:
                break
        if not nucleobase_attached:
            continue
        
        # Count hydroxyl groups (–OH) attached directly to ring atoms.
        hydroxyl_count = 0
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring:
                    continue
                if nbr.GetAtomicNum() == 8:
                    # Check if this oxygen has at least one hydrogen attached.
                    # Using GetTotalNumHs(including implicit Hs) works since we used AddHs.
                    if nbr.GetTotalNumHs() >= 1:
                        hydroxyl_count += 1
        # In many ribonucleosides the ribose ring carries at least two –OH groups
        if hydroxyl_count < 2:
            continue
        
        # Look for an exocyclic CH2OH group attached to one of the ring carbons.
        ch2oh_found = False
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            # Check only carbon atoms; they might have substituents including CH2OH.
            if atom.GetAtomicNum() != 6:
                continue
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring:
                    continue
                # Looking for a carbon substituent that could be the CH2 group
                if nbr.GetAtomicNum() == 6:
                    # A CH2 group typically has 2 bound hydrogens.
                    if nbr.GetTotalNumHs() != 2:
                        continue
                    # Now check that this carbon is attached to an OH.
                    for subnbr in nbr.GetNeighbors():
                        if subnbr.GetIdx() == atom.GetIdx():
                            continue
                        if subnbr.GetAtomicNum() == 8 and subnbr.GetTotalNumHs() >= 1:
                            ch2oh_found = True
                            break
                    if ch2oh_found:
                        break
            if ch2oh_found:
                break
        
        if not ch2oh_found:
            # If no CH2OH substituent is detected, reject this candidate.
            continue
        
        # If we reach here for this ring, then we have found a candidate ribose:
        # - five-membered ring with 4 C and 1 O,
        # - no adjacent phosphate,
        # - at least one ring carbon is bound to a nitrogen (nucleobase),
        # - at least two –OH groups on the ring,
        # - and at least one exocyclic CH2OH.
        candidate_found = True
        break

    if candidate_found:
        return True, "Contains D-ribose sugar moiety with a nucleobase attached"
    else:
        return False, "No valid D-ribose nucleoside pattern detected"

# Example testing (uncomment to run):
# test_smiles_list = [
#     "COC(=O)Cc1cn([C@@H]2O[C@H](CO)[C@@H](O)[C@H]2O)c(=O)[nH]c1=O",  # 5-methoxycarbonylmethyluridine (expected True)
#     "OC[C@H]1O[C@H]([C@H](O)[C@@H]1O)N1CCC(=O)NC1=O",                   # dihydrouridine (expected True)
#     "Nc1ncnc2n(cnc12)[C@@H]1O[C@H](CO)[C@@H](O)[C@H]1O",                   # adenosine (expected True)
#     "P1(OC2[C@](O1)([C@@H](O[C@@H]2CO)N3C=CC(=NC3=O)N)[H])(O)=O"           # example containing phosphate (expected False)
# ]
# for s in test_smiles_list:
#     result, reason = is_ribonucleoside(s)
#     print(result, reason)