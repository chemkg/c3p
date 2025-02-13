"""
Classifies: CHEBI:68452 azole
"""
"""
Classifies: Azole
Definition: Any monocyclic heteroarene consisting of a five-membered ring containing nitrogen.
Azoles can also contain one or more other non-carbon atoms, such as nitrogen, sulfur or oxygen.
"""

from rdkit import Chem

def is_azole(smiles: str):
    """
    Determines if a molecule is an azole based on its SMILES string.
    An azole is defined as any monocyclic heteroarene consisting of a five-membered ring containing nitrogen.
    The five-membered ring must be aromatic and not fused (i.e., its atoms are not shared with another ring).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule contains a qualifying azole ring, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string using RDKit
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Ensure aromaticity is perceived
    Chem.SanitizeMol(mol)

    # Get the ring information from the molecule
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    if not atom_rings:
        return False, "No rings found in the molecule"
    
    # For each ring in the molecule
    for ring in atom_rings:
        # We are only interested in rings that are exactly 5 atoms
        if len(ring) != 5:
            continue

        # Check if this ring is aromatic and non-fused:
        #   - All atoms in the five-membered ring should be aromatic.
        #   - Each atom in the ring should only belong to this ring (i.e., not be part of another ring).
        all_aromatic = True
        ring_isolated = True
        has_nitrogen = False
        
        # Count ring memberships for each atom
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if not atom.GetIsAromatic():
                all_aromatic = False
                break
            # Check if the atom is nitrogen (atomic num 7)
            if atom.GetAtomicNum() == 7:
                has_nitrogen = True
            # Count how many rings this atom participates in.
            if len(ring_info.AtomRings()[idx:idx+1]) > 1:
                # This is a bit tricky: RDKit's GetRingInfo().AtomRings() returns a tuple of rings,
                # so to check fusion, we can count occurrences of the atom index in all rings.
                count = sum(1 for r in atom_rings if idx in r)
                if count > 1:
                    ring_isolated = False
                    break
        
        if not all_aromatic:
            continue  # not an aromatic ring, skip
        if not ring_isolated:
            continue  # fused ring, skip
        if not has_nitrogen:
            continue  # does not contain nitrogen, skip

        # If we reach here, we found a valid five-membered aromatic heterocycle (azole)
        return True, "Found five-membered aromatic non-fused ring containing nitrogen (azole) in the molecule"
    
    return False, "No five-membered monocyclic aromatic ring containing nitrogen (azole) found"