"""
Classifies: CHEBI:26195 polyphenol
"""
"""
Classifies: Polyphenol
Definition: Members of the class of phenols that contain 2 or more benzene rings 
            (i.e. six-membered aromatic rings made exclusively of carbons) each of which is substituted by at least one hydroxy (-OH) group.
            
Improvements: 
 - Add explicit hydrogens so that -OH groups become detectable.
 - Only count rings that are six-membered, completely aromatic (all atoms aromatic carbons) 
   and where at least one ring atom has a neighboring oxygen that is a genuine hydroxy group 
   (i.e. oxygen bonded only to the ring atom and at least one hydrogen).
"""
from rdkit import Chem

def is_polyphenol(smiles: str):
    """
    Determines if a molecule is a polyphenol based on its SMILES string.
    
    A polyphenol is defined as a molecule that contains 2 or more benzene rings
    (six-membered aromatic rings consisting of carbon atoms) where each ring is substituted 
    by at least one hydroxy (-OH) group attached directly to one of the ring atoms.
    
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
    
    try:
        Chem.SanitizeMol(mol)
    except Exception as e:
        return False, f"Sanitization error: {e}"
    
    # Add explicit hydrogens so that hydroxy groups become explicit.
    mol = Chem.AddHs(mol)
    
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    
    # To avoid double counting, we will keep track of rings that qualify.
    phenol_ring_count = 0
    
    # Iterate over each ring of atoms in the molecule.
    for ring in atom_rings:
        # Consider only 6-membered rings.
        if len(ring) != 6:
            continue
        
        # Ensure every atom in the ring is aromatic and is a carbon.
        if not all(mol.GetAtomWithIdx(idx).GetIsAromatic() and mol.GetAtomWithIdx(idx).GetAtomicNum() == 6 for idx in ring):
            continue
        
        # For the ring to qualify, at least one of its atoms must be substituted by a hydroxy group.
        # A genuine hydroxy group here is one where an oxygen (atomic num==8) is attached to the ring atom,
        # the oxygen has only one bond (to the ring atom) and has at least one hydrogen attached.
        has_hydroxy = False
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            for neighbor in atom.GetNeighbors():
                # Skip if neighbor belongs to the current ring.
                if neighbor.GetIdx() in ring:
                    continue
                # Check if the neighbor is an oxygen.
                if neighbor.GetAtomicNum() == 8:
                    # Check that this oxygen is only bonded (degree == 1) to the ring atom:
                    # (A hydroxy oxygen should have degree 1 (to the carbon) plus its attached hydrogen(s)).
                    if neighbor.GetDegree() == 1:
                        # Also, ensure that the oxygen has at least one hydrogen (explicitly added)
                        if neighbor.GetTotalNumHs() >= 1:
                            has_hydroxy = True
                            break
            if has_hydroxy:
                break
        
        if has_hydroxy:
            phenol_ring_count += 1
            
    # There must be at least 2 benzene (phenolic) rings in the molecule for it to be classified as a polyphenol.
    if phenol_ring_count >= 2:
        return True, f"Contains {phenol_ring_count} benzene rings each substituted with a hydroxy group."
    else:
        return False, f"Contains only {phenol_ring_count} benzene ring(s) with a hydroxy substituent (need at least 2)."

# Example usage (this block can be commented out in production):
if __name__ == "__main__":
    # Example: 1,8-dihydroxy-3-methoxy-6-methylxanthone
    test_smiles = "COc1cc(O)c2c(c1)oc1cc(C)cc(O)c1c2=O"  
    result, reason = is_polyphenol(test_smiles)
    print("Is polyphenol:", result)
    print("Reason:", reason)