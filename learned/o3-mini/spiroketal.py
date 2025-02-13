"""
Classifies: CHEBI:72600 spiroketal
"""
"""
Classifies: Spiroketal 
Definition: A cyclic ketal in which the ketal carbon is the only common atom of two rings.
"""

from rdkit import Chem

def is_spiroketal(smiles: str):
    """
    Determines if a molecule is a spiroketal based on its SMILES string.
    A spiroketal is defined as a cyclic ketal where the ketal carbon is the only atom shared between two rings.
    
    Approach:
      1. Parse the SMILES string.
      2. Use RDKit to retrieve ring information.
      3. For each carbon atom present in at least two rings, check for a spiro connectivity:
         - For each pair of rings containing that atom, if the only common atom is the candidate, it is a spiro center.
      4. To qualify as a ketal center, the carbon should have at least two oxygen neighbors.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a spiroketal, False otherwise.
        str: Reason for classification.
    """
    # Parse the input SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()  # Tuple of tuples; each is a ring defined by atom indices

    # For each atom, determine which rings it is in.
    # We create a dictionary: atom index -> list of rings (as sets) that include that atom.
    atom_to_rings = {}
    for ring in atom_rings:
        ring_set = set(ring)
        for idx in ring:
            if idx not in atom_to_rings:
                atom_to_rings[idx] = []
            atom_to_rings[idx].append(ring_set)

    # Now check every carbon atom that is in at least two rings.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:
            continue  # candidate ketal center must be a carbon
        idx = atom.GetIdx()
        if idx not in atom_to_rings or len(atom_to_rings[idx]) < 2:
            continue  # need at least two rings for a spiro center

        rings_of_atom = atom_to_rings[idx]
        # For each pair of rings that contain this atom, check if intersection is exactly the spiro center.
        for i in range(len(rings_of_atom)):
            for j in range(i+1, len(rings_of_atom)):
                ring1 = rings_of_atom[i]
                ring2 = rings_of_atom[j]
                # Intersection should be exactly the candidate atom for a spiro connection.
                if ring1.intersection(ring2) == {idx}:
                    # Now verify that the candidate carbon has at least 2 oxygen neighbors 
                    # (a requirement for being the ketal center).
                    oxy_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 8]
                    if len(oxy_neighbors) >= 2:
                        return True, f"Found spiroketal center at atom index {idx} shared by two rings with at least 2 oxygen substituents."
                    else:
                        return False, f"Found spiro center at atom index {idx} but it does not have enough oxygen substituents (found {len(oxy_neighbors)})."
                    
    return False, "No spiroketal center identified (no carbon found that is the sole common atom between two rings with ketal (oxygen) substituents)."