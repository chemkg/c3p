"""
Classifies: CHEBI:18254 ribonucleoside
"""
"""
Classifies: Any nucleoside where the sugar component is D-ribose.
A ribonucleoside is defined as a nucleoside in which the sugar part is D‑ribose – a five‐membered (furanose) ring
consisting of 4 carbons and 1 oxygen, bearing hydroxyl substituents (including an exocyclic CH2OH group) and attached via a glycosidic bond to a nucleobase.
This implementation looks for:
  - A 5‑membered ring with exactly 4 carbons and 1 oxygen.
  - No adjacent phosphate group.
  - At least two –OH groups on ring atoms.
  - An exocyclic CH2OH substituent on one of the ring carbons.
  - And a ring carbon with a bond to a non‐ring nitrogen (as a rough proxy for the nucleobase attachment).
If these criteria are met the molecule is classified as a ribonucleoside.
"""
from rdkit import Chem

def is_ribonucleoside(smiles: str):
    """
    Determines if the molecule (given by its SMILES string) is a ribonucleoside,
    i.e. it contains a D‑ribose sugar moiety (a furanose ring with proper hydroxyl and CH2OH substituents)
    linked to a nucleobase (indicated by a non‐sugar nitrogen attachment).

    Args:
        smiles (str): SMILES string representing the molecule.

    Returns:
        bool: True if classified as a ribonucleoside, False otherwise.
        str: Explanation for the decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens, ensuring we count OH groups properly.
    mol = Chem.AddHs(mol)
    
    # Get all rings in the molecule.
    rings = mol.GetRingInfo().AtomRings()
    candidate_found = False
    
    for ring in rings:
        # We require a 5-membered ring.
        if len(ring) != 5:
            continue
        
        # Count atoms in the ring.
        carbon_count = 0
        oxygen_count = 0
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() == 6:
                carbon_count += 1
            elif atom.GetAtomicNum() == 8:
                oxygen_count += 1
        # A ribofuranose ring should have exactly four carbons and one oxygen.
        if carbon_count != 4 or oxygen_count != 1:
            continue
        
        # Exclude rings that are attached to a phosphate group.
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
        
        # Check for the glycosidic bond: at least one ring carbon should be bonded to a nitrogen (nucleobase proxy)
        nucleobase_attached = False
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() != 6:
                continue
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring:
                    continue
                if nbr.GetAtomicNum() == 7:  # nitrogen
                    nucleobase_attached = True
                    break
            if nucleobase_attached:
                break
        if not nucleobase_attached:
            continue
        
        # Count hydroxyl (-OH) groups directly attached to ring atoms.
        hydroxyl_count = 0
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring:
                    continue
                if nbr.GetAtomicNum() == 8:  # oxygen
                    # Check if oxygen has at least one hydrogen (indicating an OH)
                    if nbr.GetTotalNumHs() >= 1:
                        hydroxyl_count += 1
        if hydroxyl_count < 2:
            continue
        
        # Look for an exocyclic CH2OH substituent.
        ch2oh_found = False
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            # Focus on ring carbons
            if atom.GetAtomicNum() != 6:
                continue
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring:
                    continue
                if nbr.GetAtomicNum() == 6:  # candidate for CH2 group
                    # Relax: require that this carbon has at least 2 hydrogens (CH2 or even CH3 in some depictions)
                    if nbr.GetTotalNumHs() < 2:
                        continue
                    # Verify that this carbon is attached to an OH group.
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
            continue

        # If we reach here for one of the rings, we have found a candidate D-ribose moiety.
        candidate_found = True
        break

    if candidate_found:
        return True, "Contains D-ribose sugar moiety with nucleobase attachment"
    else:
        return False, "No valid D-ribose nucleoside pattern detected"

# Example tests (to run, uncomment the following lines):
# test_smiles_list = [
#     "COC(=O)Cc1cn([C@@H]2O[C@H](CO)[C@@H](O)[C@H]2O)c(=O)[nH]c1=O",  # 5-methoxycarbonylmethyluridine
#     "OC[C@H]1O[C@H]([C@H](O)[C@@H]1O)N1CCC(=O)NC1=O",                   # dihydrouridine
#     "CNC(=O)[C@H]1O[C@H]([C@H](O)[C@@H]1O)n1cnc2c(NCc3cccc(I)c3)ncnc12",    # 3-iodobenzyl-5'-N-methylcarboxamidoadenosine
#     "Nc1ncnc2n(cnc12)[C@@H]1O[C@H](CO)[C@@H](O)[C@H]1O"                    # adenosine
# ]
# for sm in test_smiles_list:
#     res, reason = is_ribonucleoside(sm)
#     print(sm, "=>", res, reason)