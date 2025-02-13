"""
Classifies: CHEBI:39362 mononitrophenol
"""
#!/usr/bin/env python3
"""
Classifies: mononitrophenol
Definition: A mononitrophenol is defined as a phenol (a benzene ring that bears at least one neutral -OH substituent)
with exactly one nitro group ([N+](=O)[O-]) directly attached to one of its aromatic carbons.
The molecule is assumed to have a molecular weight below 500 Da.
"""

from rdkit import Chem
from rdkit.Chem import Descriptors

def is_mononitrophenol(smiles: str):
    """
    Determines if a molecule is a mononitrophenol based on its SMILES string.

    The molecule must:
      - Parse correctly.
      - Have a relatively low molecular weight (< 500 Da).
      - Contain exactly one nitro group ([N+](=O)[O-]) in the molecule.
      - Have that nitro group attached directly to an aromatic carbon.
      - Possess a benzene ring (6-membered, aromatic, all-carbons) that bears at least one neutral hydroxyl substituent (-OH)
        (attached by a single bond) on one of its carbons – the same ring which bears the nitro substituent.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is classified as a mononitrophenol, False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse SMILES string into molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check molecular weight (< 500 Da)
    mol_wt = Descriptors.ExactMolWt(mol)
    if mol_wt > 500:
        return False, f"Molecular weight too high ({mol_wt:.1f} Da) for a mononitrophenol"
    
    # Define nitro group pattern.
    nitro_pattern = Chem.MolFromSmarts("[N+](=O)[O-]")
    nitro_matches = mol.GetSubstructMatches(nitro_pattern)
    if len(nitro_matches) == 0:
        return False, "No nitro group found"
    if len(nitro_matches) > 1:
        return False, f"Multiple nitro groups found ({len(nitro_matches)}), expected exactly one"
    
    # Identify the nitro group and check it is attached to an aromatic carbon.
    nitro_match = nitro_matches[0]
    # The first index in the nitro match is the nitrogen.
    nitro_N = mol.GetAtomWithIdx(nitro_match[0])
    attached_aromatic_carbons = []
    for nbr in nitro_N.GetNeighbors():
        # Skip if neighbor is part of the nitro group itself.
        if nbr.GetIdx() in nitro_match:
            continue
        # If neighbor is carbon and aromatic, register it.
        if nbr.GetAtomicNum() == 6 and nbr.GetIsAromatic():
            attached_aromatic_carbons.append(nbr.GetIdx())
    if not attached_aromatic_carbons:
        return False, "Nitro group is not attached to an aromatic carbon"
    
    # Identify candidate benzene rings.
    # A benzene ring is defined as a 6-membered ring with all atoms as carbons that are aromatic.
    ring_info = mol.GetRingInfo()
    candidate_rings = []
    for ring in ring_info.AtomRings():
        if len(ring) != 6:
            continue
        # Check that all atoms in the ring are carbon and aromatic.
        if all(mol.GetAtomWithIdx(idx).GetAtomicNum() == 6 and mol.GetAtomWithIdx(idx).GetIsAromatic() 
               for idx in ring):
            candidate_rings.append(ring)
    
    if not candidate_rings:
        return False, "No benzene ring found"
    
    # For each candidate ring, check if:
    # 1. One of the ring carbons is attached to the nitro group.
    # 2. The ring has at least one neutral hydroxyl substituent (-OH) attached by a single bond.
    valid_rings = []
    for ring in candidate_rings:
        # Check nitro attachment: at least one aromatic carbon (of the nitro group neighbors) must be in the ring.
        if not any(idx in ring for idx in attached_aromatic_carbons):
            continue
        # Check for a neutral hydroxyl group on this ring.
        hydroxyl_found = False
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            # Examine each neighbor of the ring carbon not in the ring.
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring:
                    continue
                # Check if neighbor is an oxygen.
                if nbr.GetAtomicNum() != 8:
                    continue
                # Check that the bond is single.
                bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
                if not bond or bond.GetBondTypeAsDouble() != 1.0:
                    continue
                # Check that the oxygen is neutral.
                if nbr.GetFormalCharge() != 0:
                    continue
                # Check that the oxygen has at least one hydrogen.
                if nbr.GetTotalNumHs() < 1:
                    continue
                # Found a hydroxyl substituent.
                hydroxyl_found = True
                break
            if hydroxyl_found:
                break
        if hydroxyl_found:
            valid_rings.append(ring)
    
    # Expect exactly one benzene ring that qualifies.
    if len(valid_rings) != 1:
        return False, f"Expected exactly one nitrophenol ring, found {len(valid_rings)}"
    
    return True, "Molecule is a mononitrophenol: benzene ring with one neutral hydroxyl substituent and a nitro group attached"

# Example usage:
if __name__ == '__main__':
    # Test with one of the examples: 4-nitrophenol.
    test_smiles = "Oc1ccc(cc1)[N+]([O-])=O"
    result, reason = is_mononitrophenol(test_smiles)
    print(result, reason)