"""
Classifies: CHEBI:33848 polycyclic arene
"""
"""
Classifies: polycyclic arene = polycyclic aromatic hydrocarbon
A polycyclic aromatic hydrocarbon is defined as a molecule with at least two fused (sharing at least one bond, 
i.e. sharing at least two atoms) aromatic rings.
"""

from rdkit import Chem

def is_polycyclic_arene(smiles: str):
    """
    Determines if a molecule is a polycyclic arene (polycyclic aromatic hydrocarbon) based on its SMILES string.
    
    A polycyclic aromatic hydrocarbon must contain at least two aromatic rings that are fused.
    Fused rings are defined here as aromatic rings that share at least two atoms (which indicates a shared bond).

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if molecule is a polycyclic arene, False otherwise.
        str: Reason for classification.
    """
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Ensure that aromaticity is perceived
    Chem.SanitizeMol(mol)
    
    # Get ring information (list of tuples of atom indices in each ring)
    ring_info = mol.GetRingInfo().AtomRings()
    if not ring_info or len(ring_info) < 2:
        return False, "Fewer than two rings found in molecule"
    
    # Filter for aromatic rings: ring is aromatic if every atom in the ring is aromatic.
    aromatic_rings = []
    for ring in ring_info:
        if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            aromatic_rings.append(set(ring))
    
    if len(aromatic_rings) < 2:
        return False, "Fewer than two aromatic rings detected"
    
    # Check for fused rings: two rings are fused if they share at least two atoms.
    fused_found = False
    num_aromatic = len(aromatic_rings)
    for i in range(num_aromatic):
        for j in range(i+1, num_aromatic):
            # Check intersection of atoms in the two rings
            shared_atoms = aromatic_rings[i].intersection(aromatic_rings[j])
            if len(shared_atoms) >= 2:
                fused_found = True
                break
        if fused_found:
            break
            
    if not fused_found:
        return False, "Aromatic rings are not fused (no pair shares at least two atoms)"
    
    return True, "Molecule is a polycyclic aromatic hydrocarbon with fused aromatic rings"

# Example usage:
if __name__ == "__main__":
    test_smiles = "c1ccc2c(c1)cc1ccc3cccc4ccc2c1c34"  # benzo[a]pyrene is an example
    result, reason = is_polycyclic_arene(test_smiles)
    print("Is polycyclic arene:", result)
    print("Reason:", reason)