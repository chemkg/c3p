"""
Classifies: CHEBI:26195 polyphenol
"""
"""
Classifies: Polyphenol
Definition: Members of the class of phenols that contain 2 or more benzene rings 
            each of which is substituted by at least one hydroxy group.
"""
from rdkit import Chem

def is_polyphenol(smiles: str):
    """
    Determines if a molecule is a polyphenol based on its SMILES string.
    
    A polyphenol is defined as a molecule that contains 2 or more benzene rings,
    where each benzene ring (i.e. six-membered aromatic ring) is substituted by at least one hydroxy group.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule meets the polyphenol criteria, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Ensure aromaticity and ring info are computed.
    Chem.SanitizeMol(mol)
    
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    
    phenol_ring_count = 0  # Count of benzene rings (6-membered aromatic) with at least one -OH substituent.
    # Iterate over all rings identified in the molecule.
    for ring in atom_rings:
        if len(ring) != 6:
            continue  # Only consider six-membered rings.
        # Check if every atom in the ring is aromatic.
        if not all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            continue
        
        # Now check for a hydroxy substituent attached to any atom in the ring.
        has_hydroxy = False
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            # Look at atoms connected to the ring atom that are not in the ring.
            for neighbor in atom.GetNeighbors():
                # Skip neighbors that are in the same ring.
                if neighbor.GetIdx() in ring:
                    continue
                # Check if the neighbor is oxygen and has at least one hydrogen.
                if neighbor.GetAtomicNum() == 8 and neighbor.GetTotalNumHs() >= 1:
                    has_hydroxy = True
                    break
            if has_hydroxy:
                break  # No need to check further within this ring.
        
        if has_hydroxy:
            phenol_ring_count += 1

    # Check the overall classification based on the count of phenolic benzene rings.
    if phenol_ring_count >= 2:
        return True, f"Contains {phenol_ring_count} benzene rings each substituted with a hydroxy group."
    else:
        return False, f"Contains only {phenol_ring_count} benzene ring(s) with a hydroxy substituent (need at least 2)."
    
# Example usage (You can remove or comment these lines if not needed):
if __name__ == "__main__":
    test_smiles = "COc1cc(O)c2c(c1)oc1c(O)cccc1c2=O"  # 1,8-dihydroxy-3-methoxy-6-methylxanthone
    result, reason = is_polyphenol(test_smiles)
    print("Is polyphenol:", result)
    print("Reason:", reason)