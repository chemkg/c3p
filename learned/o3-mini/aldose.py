"""
Classifies: CHEBI:15693 aldose
"""
"""
Classifies: Aldose – Aldehydic parent sugars (polyhydroxy aldehydes H[CH(OH)]nC(=O)H, n>=2) and their intramolecular hemiacetals.
Examples include beta-D-idopyranose, D-ribofuranose, L-erythrose, beta-D-glucose, etc.
"""

from rdkit import Chem

def is_aldose(smiles: str):
    """
    Determines if a molecule is an aldose based on its SMILES string.
    Aldoses are defined as polyhydroxy aldehydes (H[CH(OH)]nC(=O)H with n>=2) or their cyclic hemiacetal forms.
    Typically, aldoses have between 3 and 7 carbon atoms.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an aldose, False otherwise
        str: Reason for the classification decision
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Count carbon atoms
    carbons = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    nC = len(carbons)
    if nC < 3 or nC > 7:
        return False, f"Number of carbons is {nC}, which is outside the expected range for aldoses (3-7)"
    
    # Count hydroxyl (-OH) groups using SMARTS for oxygen with one hydrogen attached
    oh_smarts = Chem.MolFromSmarts("[OX2H]")
    oh_matches = mol.GetSubstructMatches(oh_smarts)
    nOH = len(oh_matches)
    if nOH < 2:
        return False, f"Found only {nOH} hydroxyl group(s); aldoses are polyhydroxy compounds"
    
    # Define SMARTS for an aldehyde group: a carbonyl carbon with one hydrogen [CX3H1](=O)
    aldehyde_smarts = Chem.MolFromSmarts("[CX3H1](=O)")
    has_aldehyde = mol.HasSubstructMatch(aldehyde_smarts)
    
    if has_aldehyde:
        return True, "Molecule has an open-chain aldehyde group and sufficient hydroxyl groups"
    
    # If not open-chain, try to detect a cyclic hemiacetal form.
    # Typical cyclic aldoses have 5- or 6-membered rings with a single ring oxygen.
    ring_info = mol.GetRingInfo().AtomRings()
    found_cyclic_aldose = False
    reason = ""
    
    for ring in ring_info:
        ring_size = len(ring)
        if ring_size not in (5, 6):
            continue  # Skip rings that are not furanose or pyranose sized
        # Count the number of oxygen atoms in the ring
        oxy_in_ring = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8)
        # Most sugar rings (furanose/pyranose) have exactly one oxygen.
        if oxy_in_ring != 1:
            continue
        
        # Check if at least one ring carbon bears an exocyclic hydroxyl group.
        # Loop over each atom in the ring that is carbon.
        hemiacetal_found = False
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() != 6:
                continue
            # Look at neighbors that are NOT in the ring
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring:
                    continue
                # Check if neighbor is an -OH group i.e. oxygen with one hydrogen
                if nbr.GetAtomicNum() == 8:
                    # Check if this oxygen has exactly one hydrogen neighbor
                    h_count = sum(1 for x in nbr.GetNeighbors() if x.GetAtomicNum() == 1)
                    if h_count == 1:
                        hemiacetal_found = True
                        break
            if hemiacetal_found:
                break
        
        if hemiacetal_found:
            found_cyclic_aldose = True
            reason = "Molecule is cyclic with a 5- or 6-membered ring (containing one oxygen) and an exocyclic hydroxyl indicating a hemiacetal center"
            break

    if found_cyclic_aldose:
        return True, reason
    else:
        return False, "No open-chain aldehyde group or cyclic hemiacetal sugar pattern was detected"
        
# Example test cases (uncomment to run)
# test_smiles = [
#     "OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O",  # beta-D-glucose
#     "[H]C(=O)[C@@H](O)[C@H](O)[C@H](O)CO",            # aldehydo-D-lyxose
#     "[H]C(=O)[C@H](O)CO"                              # D-glyceraldehyde
# ]
# for s in test_smiles:
#     result, message = is_aldose(s)
#     print(s, "->", result, ",", message)