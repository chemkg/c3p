"""
Classifies: CHEBI:18254 ribonucleoside
"""
"""
Classifies: Any nucleoside where the sugar component is D-ribose.
A ribonucleoside is defined as a nucleoside whose sugar is D-ribose – that is,
a five-membered (furanose) ring composed of 4 carbons and 1 oxygen, decorated with
multiple –OH groups (and a CH2OH moiety) and linked (via the anomeric carbon) to a nucleobase.
This algorithm identifies candidate sugar rings by their ring atom composition
and then verifies that a sugar-like motif is present (with several hydroxyl groups and no phosphate),
and that one of its carbons is linked to a nitrogen (from the nucleobase).
"""

from rdkit import Chem

def is_ribonucleoside(smiles: str):
    """
    Determines if a molecule is a ribonucleoside (i.e. a nucleoside where the sugar is D-ribose)
    based on its SMILES string. The algorithm looks for a five-membered ring containing exactly
    four carbon atoms and one oxygen; then it checks that several hydroxyl groups are attached
    (as expected for a ribose) and that one of the ring carbons is linked to a nitrogen (the nucleobase).
    Furthermore, if a phosphate (P) is attached to any ring atom the candidate is rejected.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule is classified as a ribonucleoside, otherwise False.
        str: Reason for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get all rings from the molecule using RDKit ring information.
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    
    candidate_found = False
    for ring in atom_rings:
        # We only consider rings of 5 atoms.
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
        # A proper ribose ring contains exactly 1 oxygen and 4 carbons.
        if oxygen_count != 1 or carbon_count != 4:
            continue

        # Check substituents on ring atoms for hydroxyl groups.
        # We count an -OH group when a neighbor is an oxygen that is bonded to a hydrogen.
        hydroxyl_count = 0
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            for nbr in atom.GetNeighbors():
                # Only consider neighbors that are outside the ring.
                if nbr.GetIdx() in ring:
                    continue
                if nbr.GetAtomicNum() == 8:
                    # Check if this oxygen has at least one hydrogen as neighbor.
                    for subnbr in nbr.GetNeighbors():
                        if subnbr.GetAtomicNum() == 1:
                            hydroxyl_count += 1
                            break

        # Require that at least three hydroxyl groups (or hydroxymethyl) are attached.
        # This is a heuristic to avoid cases where the ring could be part of another sugar.
        if hydroxyl_count < 3:
            continue

        # Exclude candidate rings that have a phosphate (P) attached.
        phosphate_attached = False
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring:
                    continue
                if nbr.GetAtomicNum() == 15:  # phosphorus
                    phosphate_attached = True
                    break
            if phosphate_attached:
                break
        if phosphate_attached:
            continue

        # Check that at least one ring carbon is attached to a nitrogen from outside the ring.
        # This nitrogen is likely part of the nucleobase.
        nucleobase_attached = False
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            # We check only carbons in the ring (the anomeric carbon is typically a carbon).
            if atom.GetAtomicNum() != 6:
                continue
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

        # If all criteria are met for this ring, we consider the molecule as a ribonucleoside.
        candidate_found = True
        break

    if candidate_found:
        return True, "Contains D-ribose sugar moiety with an attached nucleobase"
    else:
        return False, "No valid D-ribose nucleoside pattern detected"

# Example testing (you can uncomment these lines):
# test_smiles_list = [
#     "COC(=O)Cc1cn([C@@H]2O[C@H](CO)[C@@H](O)[C@H]2O)c(=O)[nH]c1=O",  # 5-methoxycarbonylmethyluridine (should be True)
#     "OC[C@H]1O[C@H]([C@H](O)[C@@H]1O)N1CCC(=O)NC1=O",                   # dihydrouridine (should be True)
#     "P1(OC2[C@](O1)([C@@H](O[C@@H]2CO)N3C=CC(=NC3=O)N)[H])(O)=O"           # false positive case (should be False)
# ]
# for s in test_smiles_list:
#     result, reason = is_ribonucleoside(s)
#     print(result, reason)