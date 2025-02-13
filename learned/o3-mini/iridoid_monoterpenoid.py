"""
Classifies: CHEBI:50563 iridoid monoterpenoid
"""
"""
Classifies: Iridoid Monoterpenoid
Definition: A monoterpenoid biosynthesized from isoprene which usually consists of a cyclopentane ring fused to 
a six-membered oxygen heterocycle, characteristic of iridoids. Cleavage of a bond in the cyclopentane may give rise to secoiridoids.
"""

from rdkit import Chem

def is_iridoid_monoterpenoid(smiles: str):
    """
    Determines if the given SMILES string corresponds to an iridoid monoterpenoid.
    The method checks for the presence of a bicyclic system with one cyclopentane ring (5-membered) 
    fused (sharing >=2 atoms) with a 6-membered ring that contains at least one oxygen atom.
    
    Args:
        smiles (str): The SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule contains the iridoid core, False otherwise.
        str: A message explaining the classification decision.
    """
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    ring_info = mol.GetRingInfo().AtomRings()  # tuple of tuples; each inner tuple contains atom indices in a ring
    if not ring_info:
        return False, "No rings detected in the molecule"

    # List candidate rings for cyclopentane (5-membered) and oxygen-containing six-membered rings
    cyclopentane_rings = []
    oxygen_six_membered_rings = []
    
    for ring in ring_info:
        if len(ring) == 5:
            cyclopentane_rings.append(ring)
        elif len(ring) == 6:
            # Check for oxygen atom in the six-membered ring
            has_oxygen = any(mol.GetAtomWithIdx(idx).GetAtomicNum() == 8 for idx in ring)
            if has_oxygen:
                oxygen_six_membered_rings.append(ring)
    
    if not cyclopentane_rings:
        return False, "No cyclopentane (5-membered) rings detected"
    if not oxygen_six_membered_rings:
        return False, "No oxygen-containing six-membered rings detected"

    # Look for a fused system: the cyclopentane and six-membered rings must share at least 2 atoms.
    for ring5 in cyclopentane_rings:
        for ring6 in oxygen_six_membered_rings:
            shared_atoms = set(ring5).intersection(set(ring6))
            if len(shared_atoms) >= 2:
                return True, "Found fused cyclopentane and oxygen heterocycle ring system characteristic of iridoid monoterpenoids"
    
    return False, "No fused bicyclic system (5-membered ring fused to a 6-membered oxygen heterocycle) found"