"""
Classifies: CHEBI:33848 polycyclic arene
"""
"""
Classifies: CHEBI:33853 polycyclic arene
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_polycyclic_arene(smiles: str):
    """
    Determines if a molecule is a polycyclic arene based on its SMILES string.
    A polycyclic arene is a polycyclic aromatic hydrocarbon consisting only of carbon atoms.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polycyclic arene, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check if all atoms are carbon (polycyclic arenes are pure hydrocarbons)
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:
            return False, "Contains non-carbon atoms"

    # Get all rings
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    
    # Check if there are at least 2 rings
    if len(rings) < 2:
        return False, "Not enough rings to be polycyclic"

    # Check if rings are fused and aromatic
    aromatic_rings = 0
    for ring in rings:
        # Check if all atoms in the ring are aromatic
        if all(mol.GetAtomWithIdx(atom).GetIsAromatic() for atom in ring):
            aromatic_rings += 1
            # Check if this ring shares at least 2 atoms with another ring
            shared_atoms = 0
            for other_ring in rings:
                if ring != other_ring:
                    shared_atoms = len(set(ring).intersection(other_ring))
                    if shared_atoms >= 2:
                        break
            else:
                return False, "Rings are not properly fused"

    # Need at least 2 aromatic rings
    if aromatic_rings < 2:
        return False, f"Only {aromatic_rings} aromatic ring(s) found, need at least 2"

    return True, "Contains multiple fused aromatic rings consisting only of carbon atoms"