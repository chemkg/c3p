"""
Classifies: CHEBI:75874 para-terphenyl
"""
"""
Classifies: Para-terphenyl class.
Definition: A ring assembly based on a 1,4-diphenylbenzene skeleton and its substituted derivatives.
That is, a central benzene ring substituted (in non-fused fashion) with two benzene rings at para positions.
"""
from rdkit import Chem

def is_para_terphenyl(smiles: str):
    """
    Determines if a molecule belongs to the para-terphenyl class based on its SMILES string.
    The criteria is that one of the benzene rings (six-membered aromatic rings) in the molecule
    must be substituted at two para positions with two distinct benzene rings
    (with which it shares exactly one atom via a single bond). This corresponds to a 1,4-diphenylbenzene core.
    
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
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    
    # Collect candidate benzene (six-membered aromatic) rings.
    benzene_rings = []
    for ring in atom_rings:
        if len(ring) == 6 and all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            benzene_rings.append(ring)
    
    if len(benzene_rings) < 3:
        return False, "Fewer than three benzene rings detected. A para-terphenyl core requires a central benzene and two peripheral benzene rings."
    
    # For each candidate benzene ring, treat it as the potential central ring.
    # We want to see if it has two attachments (via a single non-fused bond) to two other benzene rings,
    # and that these attachments occur at para positions (opposite vertices of the ring).
    for central_ring in benzene_rings:
        # For clarity, we assume the order provided in the ring (central_ring) is the ring ordering.
        connections = []  # list of tuples: (position_index, peripheral_ring index)
        # Loop over each atom (by order) in the central ring
        for pos, atom_idx in enumerate(central_ring):
            atom = mol.GetAtomWithIdx(atom_idx)
            # Look at each neighbor that is NOT in the central ring.
            for nb in atom.GetNeighbors():
                nb_idx = nb.GetIdx()
                if nb_idx in central_ring:
                    continue  # skip neighbors that are in the central ring itself
                # For every benzene ring candidate (peripheral candidate), check if nb is one of its atoms.
                for per_idx, per_ring in enumerate(benzene_rings):
                    # We do not want to consider the central ring itself as peripheral.
                    if per_ring == central_ring:
                        continue
                    # Check if the peripheral ring and the central ring share exactly one atom.
                    # (This avoids fused systems.)
                    common = set(central_ring).intersection(set(per_ring))
                    if len(common) != 1:
                        continue
                    # Check that the bond from the central atom to the neighbor is a single (or aromatic) bond.
                    for bond in mol.GetAtomWithIdx(atom_idx).GetBonds():
                        if bond.GetOtherAtom(atom).GetIdx() == nb_idx:
                            if bond.GetBondType() == Chem.rdchem.BondType.SINGLE or bond.GetIsAromatic():
                                # We found an attachment from the central ring (at position pos) to a peripheral benzene ring.
                                # Record the peripheral ring id (we use its index in benzene_rings to differentiate).
                                connections.append((pos, per_idx))
                            break
        # Now, in our central ring candidate, we require at least two distinct attachments.
        # Moreover the positions of attachment must be para relative; for a benzene ring, para positions 
        # are separated by 3 edges (in a 6-membered ring, the minimum of |i-j| and 6-|i-j| should be 3).
        if len(connections) < 2:
            continue
        # Check each distinct pair among the connections.
        for i in range(len(connections)):
            for j in range(i+1, len(connections)):
                pos_i, per_id_i = connections[i]
                pos_j, per_id_j = connections[j]
                # They must be attached to two distinct benzene rings.
                if per_id_i == per_id_j:
                    continue
                # Compute circular distance on a six-membered ring.
                diff = abs(pos_i - pos_j)
                circ_dist = min(diff, 6 - diff)
                if circ_dist == 3:
                    return True, "Para-terphenyl core found: central benzene ring with two para-connected benzene rings."
                    
    return False, "Para-terphenyl core not found in the molecule."

# Debug/test section (can be removed or commented out in production)
if __name__ == "__main__":
    # Some test SMILES strings:
    test_smiles_true = "c1ccc(cc1)-c2ccc(cc2)-c3ccccc3"  # Basic 1,4-diphenylbenzene example.
    test_smiles_false = "c1ccccc1"  # Just benzene.
    res_t, reason_t = is_para_terphenyl(test_smiles_true)
    res_f, reason_f = is_para_terphenyl(test_smiles_false)
    print("Test true:", res_t, reason_t)
    print("Test false:", res_f, reason_f)