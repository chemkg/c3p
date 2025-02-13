"""
Classifies: CHEBI:39362 mononitrophenol
"""
#!/usr/bin/env python3
"""
Classifies: mononitrophenol
Definition: A nitrophenol is defined as a phenol (an aromatic ring bearing a hydroxyl group in its neutral form)
that carries exactly one nitro substituent ([N+](=O)[O-]) directly attached to an aromatic carbon.
In addition, the molecule should be “small” (as most mononitrophenols are) so that large, multi‐ring compounds
with a nitrophenol substructure are not mis‐classified.
"""

from rdkit import Chem
from rdkit.Chem import Descriptors

def is_mononitrophenol(smiles: str):
    """
    Determines if a molecule is a mononitrophenol based on its SMILES string.
    
    The molecule must:
      - Parse correctly
      - Have a relatively low molecular weight (< ~500 Da)
      - Contain exactly one nitro group ([N+](=O)[O-])
      - Have that nitro group attached directly to an aromatic carbon in a six-membered (benzene-like) ring
      - Have at least one aromatic hydroxyl substituent on the same ring.
        Here we require that an oxygen is attached by a single bond and is in the neutral state;
        this helps avoid matching deprotonated phenolates (which have been troublesome).
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is classified as a mononitrophenol, False otherwise.
        str: Explanation of the classification decision.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Reject molecules that are too big; mononitrophenols are usually small.
    mol_wt = Descriptors.ExactMolWt(mol)
    if mol_wt > 500:
        return False, f"Molecular weight too high ({mol_wt:.1f} Da) for a mononitrophenol"

    # Define a pattern for the nitro group.
    nitro_pattern = Chem.MolFromSmarts("[N+](=O)[O-]")
    nitro_matches = mol.GetSubstructMatches(nitro_pattern)
    if len(nitro_matches) == 0:
        return False, "No nitro group found"
    if len(nitro_matches) > 1:
        return False, f"Multiple nitro groups found ({len(nitro_matches)}), expected exactly one"

    # Find the aromatic carbon(s) that attach(es) the nitro group.
    nitro_match = nitro_matches[0]
    nitro_N = mol.GetAtomWithIdx(nitro_match[0])
    nitro_attached_carbons = []
    for neighbor in nitro_N.GetNeighbors():
        # Skip atoms that are part of the nitro group.
        if neighbor.GetIdx() in nitro_match:
            continue
        # Check for an aromatic carbon.
        if neighbor.GetAtomicNum() == 6 and neighbor.GetIsAromatic():
            nitro_attached_carbons.append(neighbor.GetIdx())
    if not nitro_attached_carbons:
        return False, "Nitro group is not attached to an aromatic carbon"
    
    # Look for candidate rings: six-membered aromatic rings that are benzene-like.
    ring_info = mol.GetRingInfo()
    candidate_rings = []
    for ring in ring_info.AtomRings():
        if len(ring) != 6:
            continue
        # Check that every atom in the ring is aromatic and that at least 4 are carbons.
        aromatic_atoms = 0
        carbon_atoms = 0
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetIsAromatic():
                aromatic_atoms += 1
            if atom.GetAtomicNum() == 6:
                carbon_atoms += 1
        if aromatic_atoms == 6 and carbon_atoms >= 4:
            candidate_rings.append(ring)

    # From the candidate rings, select those that
    # (a) contain at least one of the nitro-attached aromatic carbons, and
    # (b) have at least one “phenol” substituent (an oxygen attached by a single bond in a neutral state).
    valid_rings = []
    for ring in candidate_rings:
        if not any(idx in ring for idx in nitro_attached_carbons):
            continue  # Nitro group is not on this ring.
        hydroxyl_found = False
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            # Look at neighbors that are not part of the ring.
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring:
                    continue
                if nbr.GetAtomicNum() == 8:
                    # Check that the bond is a single bond.
                    bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
                    if bond is None or bond.GetBondTypeAsDouble() != 1.0:
                        continue
                    # We expect a genuine hydroxyl: the oxygen should have only one heavy atom neighbor—
                    # so that it is not part of a carbonyl/carboxylate—and carry no formal charge.
                    heavy_neighbors = [a for a in nbr.GetNeighbors() if a.GetAtomicNum() != 1]
                    if len(heavy_neighbors) != 1:
                        continue
                    if nbr.GetFormalCharge() != 0:
                        continue
                    hydroxyl_found = True
                    break
            if hydroxyl_found:
                break
        if hydroxyl_found:
            valid_rings.append(ring)

    if len(valid_rings) != 1:
        return (False, f"Expected exactly one nitrophenol ring, found {len(valid_rings)}")
    
    return (True, "Molecule is a mononitrophenol: contains a phenol ring with a single nitro substituent attached to an aromatic carbon")

# Example usage:
if __name__ == '__main__':
    # Test with one of the examples: 4-nitrophenol
    test_smiles = "Oc1ccc(cc1)[N+]([O-])=O"
    result, reason = is_mononitrophenol(test_smiles)
    print(result, reason)