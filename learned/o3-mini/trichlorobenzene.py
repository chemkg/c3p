"""
Classifies: CHEBI:27096 trichlorobenzene
"""
"""
Classifies: trichlorobenzene
Definition: Any chlorobenzene (benzene ring) carrying exactly three chloro substituents at unspecified positions.
We now require that the benzene ring is isolated (i.e. not fused with additional rings via two or more shared atoms)
so as to avoid false positives from polycyclic aromatic compounds.
"""

from rdkit import Chem

def is_trichlorobenzene(smiles: str):
    """
    Determines if a molecule qualifies as a trichlorobenzene based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule contains an isolated benzene ring with exactly three chlorine substituents.
        str: Explanation of the classification decision.
    """
    # Parse the SMILES string into an RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Sanitize to ensure aromaticity is correctly perceived
    try:
        Chem.SanitizeMol(mol)
    except Exception as e:
        return False, f"Sanitization failed: {e}"
    
    # Get ring information from the molecule.
    # This returns a tuple of tuples where each inner tuple contains the indices of atoms in that ring.
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    
    # Loop over each candidate ring with 6 atoms (potential benzene ring)
    for ring in atom_rings:
        # Only consider six-membered rings
        if len(ring) != 6:
            continue
        
        # Verify that every atom in the ring is an aromatic carbon
        is_benzene = True
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() != 6 or not atom.GetIsAromatic():
                is_benzene = False
                break
        if not is_benzene:
            continue
        
        # Check whether the candidate ring is fused.
        # We define a ring as fused if it shares two or more atoms with any other ring in the molecule.
        fused = False
        for other_ring in atom_rings:
            if other_ring == ring:
                continue
            # Calculate the intersection of the two rings
            common_atoms = set(ring).intersection(other_ring)
            if len(common_atoms) >= 2:
                fused = True
                break
        # If the ring is fused, then skip it (we want an isolated benzene ring)
        if fused:
            continue

        # Count chlorine substituents attached directly to ring atoms (neighbors not in the ring)
        cl_count = 0
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            for nbr in atom.GetNeighbors():
                # Ignore neighbors that are part of the candidate ring
                if nbr.GetIdx() in ring:
                    continue
                # Count if the neighbor is a chlorine (atomic number 17)
                if nbr.GetAtomicNum() == 17:
                    cl_count += 1
        
        # If exactly 3 chlorine substituents are attached, we classify the molecule as a trichlorobenzene.
        if cl_count == 3:
            return True, "Found an isolated benzene ring with exactly three chloro substituents"
    
    return False, "No isolated benzene ring with exactly three chloro substituents found"