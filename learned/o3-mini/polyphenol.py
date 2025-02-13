"""
Classifies: CHEBI:26195 polyphenol
"""
"""
Classifies: Polyphenol
Definition: A polyphenol is a phenolic compound with 2 or more benzene rings—
(each benzene ring being a six‐membered aromatic ring made exclusively of carbon atoms)
—and each such ring has at least one hydroxy (-OH) group directly attached.
"""
from rdkit import Chem

def is_polyphenol(smiles: str):
    """
    Determines if a molecule is a polyphenol based on its SMILES string.
    
    A polyphenol is defined as a molecule that contains 2 or more benzene rings—
    six-membered aromatic rings composed solely of carbon atoms—
    where each such ring has at least one directly attached hydroxy (-OH) group.
    
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
    
    # Retrieve ring information from the molecule.
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    
    qualifying_rings = 0  # Count for benzene rings with an -OH substituent.
    
    for ring in atom_rings:
        # Consider only rings of size 6.
        if len(ring) != 6:
            continue
        
        # Check that all atoms in the ring are carbon and aromatic.
        if not all(mol.GetAtomWithIdx(idx).GetAtomicNum() == 6 and 
                   mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            continue
        
        # Now check if at least one ring carbon has a hydroxy substituent.
        has_hydroxy = False
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            # Look at neighbors that are not part of the ring.
            for neighbor in atom.GetNeighbors():
                if neighbor.GetIdx() in ring:
                    continue
                # A hydroxy group is represented by an oxygen with one or more hydrogens.
                if neighbor.GetAtomicNum() == 8 and neighbor.GetTotalNumHs() >= 1:
                    has_hydroxy = True
                    break
            if has_hydroxy:
                break
        
        if has_hydroxy:
            qualifying_rings += 1

    if qualifying_rings >= 2:
        return True, f"Contains {qualifying_rings} benzene ring(s) each substituted with a hydroxy group."
    else:
        return False, f"Contains only {qualifying_rings} benzene ring(s) with a hydroxy substituent (need at least 2)."

# Example usage:
if __name__ == "__main__":
    # Test one of the example SMILES:
    test_smiles = "COc1cc(O)c2c(c1)oc1cc(C)cc(O)c1c2=O"  # 1,8-dihydroxy-3-methoxy-6-methylxanthone
    result, reason = is_polyphenol(test_smiles)
    print("Is polyphenol:", result)
    print("Reason:", reason)