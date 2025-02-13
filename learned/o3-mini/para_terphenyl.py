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
    A para-terphenyl is defined as a ring assembly with a 1,4-diphenylbenzene skeleton 
    (i.e. a central benzene ring substituted at opposite (para) positions with two other benzene rings)
    possibly with additional substituents.
    
    This implementation first collects all six-membered aromatic rings. Then, for each pair of candidate benzene rings,
    it determines if they are connected via a single (non-fused) bond (i.e. they share exactly one atom).
    Finally, for each benzene ring we check if there are two distinct external connections at positions that are para 
    relative (separated by three atoms along the ring order).
    
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
    atom_rings = ri.AtomRings()
    
    # Collect candidate benzene rings: six atoms, all aromatic.
    benzene_rings = []  # list of tuples of atom indices
    for ring in atom_rings:
        if len(ring) == 6 and all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            benzene_rings.append(ring)
    
    if len(benzene_rings) < 3:
        return False, "Less than three aromatic six-membered rings detected. p-Terphenyl requires a 1,4-diphenylbenzene core."
    
    # For each candidate benzene ring, we will record external connections.
    # A connection is defined when an atom in one benzene ring is bonded via a single bond to an atom in a different benzene ring,
    # and the two rings share exactly one atom overall (so they are not fused).
    # We'll store for each candidate central ring a list of tuples:
    #   (pos, other_ring_index, atom idx in other ring)
    # where 'pos' is the index (position) of the atom within the ring (according to the ring ordering) that connects externally.
    connections = {i: [] for i in range(len(benzene_rings))}
    
    # Loop over each pair of candidate rings. To avoid double‐counting, we compare each pair only once.
    for i in range(len(benzene_rings)):
        ring_i = benzene_rings[i]
        set_i = set(ring_i)
        for j in range(len(benzene_rings)):
            if i == j:
                continue
            ring_j = benzene_rings[j]
            set_j = set(ring_j)
            # We want to check if ring_i and ring_j are connected via exactly one common atom.
            common = set_i.intersection(set_j)
            if len(common) != 1:
                continue  # skip if fused (more than one) or not connected at all
            common_atom = list(common)[0]
            # Now find the bond that links the two rings.
            # In ring_i, the connection atom (common_atom) might be bonded to an atom in ring_j;
            # we verify that at least one bond from common_atom connects to an atom that is in ring_j but not in ring_i except the common atom itself.
            atom_i = mol.GetAtomWithIdx(common_atom)
            found = False
            for bond in atom_i.GetBonds():
                # Look for a bond that connects an atom inside ring_i to an atom inside ring_j.
                a1 = bond.GetBeginAtom()
                a2 = bond.GetEndAtom()
                idx1, idx2 = a1.GetIdx(), a2.GetIdx()
                if (idx1 in set_i and idx2 in set_j and idx1 not in set_j) or (idx2 in set_i and idx1 in set_j and idx2 not in set_i):
                    # Verify that the bond is a single bond (the aromatic bond type in RDKit is often SINGLE with an aromatic flag).
                    if bond.GetBondType() == Chem.rdchem.BondType.SINGLE or bond.GetIsAromatic():
                        found = True
                        break
            # If found, record the connection for ring_i.
            if found:
                # In ring_i, we want to record which position corresponds to the connection.
                # The connection is at the common atom.
                try:
                    pos = ring_i.index(common_atom)
                    # Record this connection if not already recorded from ring_i to ring_j.
                    if not any(x[1] == j for x in connections[i]):
                        connections[i].append((pos, j, common_atom))
                    # Also record the reverse connection for ring_j.
                    pos_j = ring_j.index(common_atom)
                    if not any(x[1] == i for x in connections[j]):
                        connections[j].append((pos_j, i, common_atom))
                except Exception:
                    pass

    # Now search for a candidate central ring that has connections to at least 2 distinct other rings
    # and where the connection positions are in a para (i.e. 1,4-) relationship.
    for i, conn_list in connections.items():
        if len(conn_list) < 2:
            continue
        # We now examine connection positions on ring i.
        # Sort the positions (the order in ring_i as given by the ring ordering)
        positions = sorted([pos for pos, other_idx, atm in conn_list])
        # Check if any two positions are separated by 3 (or equivalently, 6-3 = 3) modulo 6.
        for a in range(len(positions)):
            for b in range(a+1, len(positions)):
                diff = abs(positions[b] - positions[a])
                circ_dist = min(diff, 6 - diff)
                if circ_dist == 3:
                    # Also confirm that the two connections come from distinct outer rings.
                    outer1 = conn_list[a][1]
                    outer2 = conn_list[b][1]
                    if outer1 != outer2:
                        return True, "Para-terphenyl core found: central benzene ring with two para-connected benzene rings."
    return False, "Para-terphenyl core not found in the molecule."

# Debug/test section (can be removed or commented out in production)
if __name__ == "__main__":
    # Some test SMILES (one from true positives and one from a false positive shown above)
    test_smiles_true = "COc1cc(-c2ccccc2)c(OC)c(O)c1-c1ccc(O)c(O)c1"  # 4''-dehydro-3-hydroxyterphenyllin
    test_smiles_false = "COC1=CC(=CC(=C1)C2=CC3=C(C=C(C=C3)OCC(=O)N)OC2=O)OC"  # a false positive example
    res_true, reason_true = is_para_terphenyl(test_smiles_true)
    res_false, reason_false = is_para_terphenyl(test_smiles_false)
    print("True test:", res_true, reason_true)
    print("False test:", res_false, reason_false)