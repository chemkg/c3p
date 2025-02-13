"""
Classifies: CHEBI:27096 trichlorobenzene
"""
"""
Classifies: trichlorobenzene
Definition: Any chlorobenzene (benzene ring) carrying three chloro substituents at unspecified positions.
The algorithm searches for any six-membered aromatic ring that is composed of carbon atoms and counts the
number of chlorine atoms directly bound (as substituents) to the ring. If at least one such ring has exactly
three chlorine substituents, the molecule is classified as a trichlorobenzene.
"""

from rdkit import Chem

def is_trichlorobenzene(smiles: str):
    """
    Determines if a molecule is a trichlorobenzene based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a trichlorobenzene, False otherwise.
        str: Explanation/reason for the classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Ensure aromaticity is properly perceived
    Chem.SanitizeMol(mol)
    
    # Retrieve ring information from the molecule
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()  # tuple of atom index tuples
    
    # Loop over each ring in the molecule
    for ring in atom_rings:
        # Only consider six-membered rings
        if len(ring) != 6:
            continue

        # Verify that every atom in the ring is aromatic carbon
        is_aromatic_benzene = True
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            # Check if the atom is carbon and aromatic
            if atom.GetAtomicNum() != 6 or not atom.GetIsAromatic():
                is_aromatic_benzene = False
                break
        if not is_aromatic_benzene:
            continue

        # Count chlorine substituents on the ring:
        # For each atom in the ring, consider its neighbors that are not part of the ring.
        cl_count = 0
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            for nbr in atom.GetNeighbors():
                # Ignore neighbors that are also in the ring
                if nbr.GetIdx() in ring:
                    continue
                # Count if neighbor is chlorine (atomic number 17)
                if nbr.GetAtomicNum() == 17:
                    cl_count += 1
        # Check if exactly three Cl substituents are attached to the ring
        if cl_count == 3:
            return True, "Found a benzene ring with exactly three chloro substituents"
    
    return False, "No benzene ring with exactly three chloro substituents found"