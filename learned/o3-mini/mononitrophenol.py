"""
Classifies: CHEBI:39362 mononitrophenol
"""
#!/usr/bin/env python3
"""
Classifies: mononitrophenol
Definition: A nitrophenol is defined as a phenol (an aromatic benzene ring bearing at least one neutral hydroxyl substituent)
that carries exactly one nitro substituent ([N+](=O)[O-]) directly attached to one of its aromatic carbons.
In addition, the molecule is assumed to be “small” (typically molecular weight <500 Da) so that large, multi‐ring compounds
are not mis‐classified.
"""

from rdkit import Chem
from rdkit.Chem import Descriptors

def is_mononitrophenol(smiles: str):
    """
    Determines if a molecule is a mononitrophenol based on its SMILES string.

    The molecule must:
      - Parse correctly.
      - Have a relatively low molecular weight (< 500 Da).
      - Contain exactly one nitro group ([N+](=O)[O-]) in the whole molecule.
      - Have that nitro group attached directly to an aromatic carbon.
      - Possess a benzene (six-membered, all-carbon aromatic) ring that bears at least one neutral hydroxyl substituent (–OH)
        (with the O attached by a single bond and having no formal charge), with the nitro group on that ring.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is classified as a mononitrophenol, False otherwise.
        str: Explanation of the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check molecular weight constraint.
    mol_wt = Descriptors.ExactMolWt(mol)
    if mol_wt > 500:
        return False, f"Molecular weight too high ({mol_wt:.1f} Da) for a mononitrophenol"
        
    # Define nitro group SMARTS.
    nitro_pattern = Chem.MolFromSmarts("[N+](=O)[O-]")
    nitro_matches = mol.GetSubstructMatches(nitro_pattern)
    if len(nitro_matches) == 0:
        return False, "No nitro group found"
    if len(nitro_matches) > 1:
        return False, f"Multiple nitro groups found ({len(nitro_matches)}), expected exactly one"
    
    # From the unique nitro match, find the nitrogen atom and then the aromatic carbon(s) it is attached to.
    nitro_match = nitro_matches[0]
    nitro_N = mol.GetAtomWithIdx(nitro_match[0])
    attached_aromatic_carbons = []
    for nbr in nitro_N.GetNeighbors():
        # Skip if neighbor is in the nitro group.
        if nbr.GetIdx() in nitro_match:
            continue
        if nbr.GetAtomicNum() == 6 and nbr.GetIsAromatic():
            attached_aromatic_carbons.append(nbr.GetIdx())
    if not attached_aromatic_carbons:
        return False, "Nitro group is not attached to an aromatic carbon"
    
    # Identify candidate benzene rings:
    # A benzene ring is defined as a 6-membered ring whose atoms are all carbon and aromatic.
    ring_info = mol.GetRingInfo()
    candidate_rings = []
    for ring in ring_info.AtomRings():
        if len(ring) != 6:
            continue
        # Check that every atom in the ring is a carbon and aromatic.
        if all(mol.GetAtomWithIdx(idx).GetAtomicNum() == 6 and mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            candidate_rings.append(ring)
    
    if not candidate_rings:
        return False, "No benzene ring found"
    
    # For identifying a neutral hydroxyl group, we use a SMARTS for -OH where oxygen is sp3 with one bond and no formal charge.
    hydroxyl_smarts = "[OX2H]"
    hydroxyl_pattern = Chem.MolFromSmarts(hydroxyl_smarts)
    
    valid_rings = []
    for ring in candidate_rings:
        # Check that one of the aromatic carbons on the ring is the attachment point of the nitro.
        if not any(idx in ring for idx in attached_aromatic_carbons):
            continue  # Nitro group not on this ring.
        # Check for at least one neutral -OH group directly attached to a carbon in the ring.
        hydroxyl_found = False
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            # For each neighbor not in the ring, look for a neutral hydroxyl group.
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring:
                    continue
                # Check if the neighbor atom matches the hydroxy pattern
                # Also verify the bond is single.
                bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
                if bond is None or bond.GetBondTypeAsDouble() != 1.0:
                    continue
                if nbr.HasSubstructMatch(hydroxyl_pattern):
                    hydroxyl_found = True
                    break
            if hydroxyl_found:
                break
        if hydroxyl_found:
            valid_rings.append(ring)
    
    # We expect exactly one benzene ring to qualify as our nitrophenol ring.
    if len(valid_rings) != 1:
        return False, f"Expected exactly one nitrophenol ring, found {len(valid_rings)}"
    
    return True, "Molecule is a mononitrophenol: benzene ring with one neutral hydroxyl substituent and one nitro group attached"

# Example usage:
if __name__ == '__main__':
    # Try one of the examples: 4-nitrophenol.
    test_smiles = "Oc1ccc(cc1)[N+]([O-])=O"
    result, reason = is_mononitrophenol(test_smiles)
    print(result, reason)