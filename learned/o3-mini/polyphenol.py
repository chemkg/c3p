"""
Classifies: CHEBI:26195 polyphenol
"""
"""
Classifies: Polyphenol
Definition: A polyphenol is a phenolic compound with 2 or more benzene rings—
(each benzene ring being a six-membered aromatic ring made exclusively of carbons)
–and each such ring has at least one hydroxy (-OH) group directly attached.
"""
from rdkit import Chem

def is_polyphenol(smiles: str):
    """
    Determines if a molecule is a polyphenol based on its SMILES string.

    A polyphenol is defined as a molecule that contains 2 or more benzene rings
    (six-membered aromatic rings consisting solely of carbon atoms) where each ring
    has at least one hydroxy (-OH) group directly attached.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule meets the polyphenol criteria, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    try:
        Chem.SanitizeMol(mol)
    except Exception as e:
        return False, f"Sanitization error: {e}"
    
    # Add explicit hydrogens so that -OH groups are explicit.
    mol = Chem.AddHs(mol)
    
    # Get the molecule's ring information.
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    
    # Count rings that qualify as benzene rings with at least one hydroxy substituent.
    qualifying_rings = 0
    for ring in atom_rings:
        # Consider only 6-membered rings.
        if len(ring) != 6:
            continue

        # Ensure all atoms in the ring are aromatic carbons (atomic num 6 and aromatic).
        if not all(mol.GetAtomWithIdx(idx).GetAtomicNum() == 6 and mol.GetAtomWithIdx(idx).GetIsAromatic() 
                   for idx in ring):
            continue
        
        # Check if this ring has at least one directly attached -OH group.
        # We iterate over the ring atoms and look at neighbors not in the ring.
        has_hydroxy = False
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            for neighbor in atom.GetNeighbors():
                # If the neighbor is part of the ring, ignore it.
                if neighbor.GetIdx() in ring:
                    continue

                # Check if the neighbor is an oxygen.
                if neighbor.GetAtomicNum() == 8:
                    # Count heavy neighbors (neighbors that are not hydrogens).
                    heavy_neighbors = sum(1 for n in neighbor.GetNeighbors() if n.GetAtomicNum() != 1)
                    # A genuine hydroxy oxygen is attached to one heavy atom (the ring carbon)
                    # and should have at least one hydrogen.
                    if heavy_neighbors == 1 and neighbor.GetTotalNumHs() >= 1:
                        has_hydroxy = True
                        break
            if has_hydroxy:
                break
        
        if has_hydroxy:
            qualifying_rings += 1
            
    if qualifying_rings >= 2:
        return True, f"Contains {qualifying_rings} benzene rings each substituted with a hydroxy group."
    else:
        return False, f"Contains only {qualifying_rings} benzene ring(s) with a hydroxy substituent (need at least 2)."

# Example usage:
if __name__ == "__main__":
    # Test with one of the example SMILES:
    test_smiles = "COc1cc(O)c2c(c1)oc1cc(C)cc(O)c1c2=O"  # 1,8-dihydroxy-3-methoxy-6-methylxanthone
    result, reason = is_polyphenol(test_smiles)
    print("Is polyphenol:", result)
    print("Reason:", reason)