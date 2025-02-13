"""
Classifies: CHEBI:75874 para-terphenyl
"""
"""
Classifies: A ring assembly based on a 1,4-diphenylbenzene skeleton and its substituted derivatives thereof – para-terphenyl.
"""
from rdkit import Chem

def is_para_terphenyl(smiles: str):
    """
    Determines if a molecule belongs to the para-terphenyl class based on its SMILES string.
    A para-terphenyl is defined as a ring assembly with a 1,4-diphenylbenzene skeleton (i.e. a central benzene ring 
    substituted at opposite (para) positions with two other benzene rings), possibly with additional substituents.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as para-terphenyl, False otherwise.
        str: Reason for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Retrieve ring information.
    ri = mol.GetRingInfo()
    # Collect rings that are six-membered and aromatic (benzene rings).
    aromatic_rings = []
    for ring in ri.AtomRings():
        if len(ring) == 6 and all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            aromatic_rings.append(ring)
    
    # At least three benzene rings are expected for a para-terphenyl skeleton.
    if len(aromatic_rings) < 3:
        return False, "Less than three aromatic six-membered rings detected. p-Terphenyl requires a 1,4-diphenylbenzene core."
    
    # We now check for a central benzene ring that has two distinct connections (bonds) to other benzene rings.
    # These connections must be at para positions (i.e. separated by 3 atoms on a six-membered ring).
    for ring in aromatic_rings:
        # List to collect positions (using ring ordering) where an atom in the candidate ring connects to another benzene ring.
        connection_positions = []
        for atom_idx in ring:
            atom = mol.GetAtomWithIdx(atom_idx)
            # Check each bond from this atom.
            for bond in atom.GetBonds():
                nbr = bond.GetOtherAtom(atom)
                # If the neighbor is in the same ring, skip.
                if nbr.GetIdx() in ring:
                    continue
                # Check if neighbor is part of any other aromatic six-membered ring.
                for other_ring in aromatic_rings:
                    if other_ring == ring:
                        continue
                    if nbr.GetIdx() in other_ring:
                        # Record the index of the atom in the ring.
                        # Note: The order provided by AtomRings() is taken as the ring ordering.
                        connection_positions.append(ring.index(atom_idx))
                        break  # Found a valid connection; move to next bond.
        # Remove duplicates (if an atom has multiple bonds to atoms in another benzene ring, count it only once)
        connection_positions = list(set(connection_positions))
        
        # For a proper 1,4-disubstitution pattern, we expect at least one pair of connections whose positions
        # differ by 3 (or 6-3 = 3) along the ring.
        if len(connection_positions) >= 2:
            connection_positions.sort()
            found_para = False
            for i in range(len(connection_positions)):
                for j in range(i+1, len(connection_positions)):
                    diff = abs(connection_positions[i] - connection_positions[j])
                    # For a cyclic ring of six, the minimum distance (accounting for wrap-around) is:
                    circ_distance = min(diff, 6 - diff)
                    if circ_distance == 3:
                        found_para = True
                        break
                if found_para:
                    break
            if found_para:
                return True, "Para-terphenyl core found: central benzene ring with two para-connected benzene rings."
    
    return False, "Para-terphenyl core not found in the molecule."