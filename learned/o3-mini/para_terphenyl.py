"""
Classifies: CHEBI:75874 para-terphenyl
"""
"""
Classifies: Para-terphenyl class.
Definition: A ring assembly based on a 1,4-diphenylbenzene skeleton and its substituted derivatives.
That is, a central benzene ring substituted (non-fused) at two para positions with two benzene rings.
Note: This implementation relaxes the connection criterion to allow one or two linker atoms
      between the central ring and a peripheral benzene.
"""

from rdkit import Chem

def is_para_terphenyl(smiles: str):
    """
    Determines if a molecule belongs to the para-terphenyl class based on its SMILES string.
    Our strategy is:
      1. Identify all six-membered aromatic (benzene) rings in the molecule.
      2. For each candidate central benzene ring, for each atom in that ring, look for a path 
         (of length up to 2 bonds) from the ring atom to any other benzene ring that is not fused
         to the central ring. Record the attachment positions together with which peripheral ring was reached.
      3. For the central ring, check if there are two distinct attachments that occur at para positions
         (i.e. positions separated by 3 in the cyclic order).
         
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as para-terphenyl, False otherwise.
        str: Reason for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get aromatic ring information (rings with 6 atoms that are all aromatic are considered benzene rings)
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    benzene_rings = []
    for ring in atom_rings:
        if len(ring) == 6 and all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            benzene_rings.append(ring)
    
    if len(benzene_rings) < 3:
        return False, "Fewer than three benzene rings detected. A para-terphenyl core requires a central benzene and two peripheral benzene rings."
    
    # Helper: given a starting atom, search up to max_depth bonds (excluding atoms belonging to central ring)
    # to see if we can reach an atom that is a member of one of the peripheral benzene rings.
    def dfs_find_peripheral(start_atom_idx, excluded, max_depth=2):
        visited = set()
        stack = [(start_atom_idx, 0)]
        while stack:
            current_idx, depth = stack.pop()
            if current_idx in visited:
                continue
            visited.add(current_idx)
            atom = mol.GetAtomWithIdx(current_idx)
            # If the current atom belongs to any benzene ring not in excluded,
            # then return that benzene ring id (index in benzene_rings list).
            for ring_id, ring in enumerate(benzene_rings):
                # Only consider rings that are different from the central ring (which we pass in "excluded")
                if ring_id in excluded:
                    continue
                if current_idx in ring:
                    return ring_id
            # If depth limit not reached then search neighbors
            if depth < max_depth:
                for nb in atom.GetNeighbors():
                    nb_idx = nb.GetIdx()
                    # Do not go back into the central ring (we use display of the central ring separately)
                    if nb_idx in excluded:
                        continue
                    stack.append((nb_idx, depth + 1))
        return None  # not found
    
    # Now, iterate over each candidate benzene ring as potential central ring.
    # We need to check if a given benzene ring (with 6 atoms) is substituted at two para positions.
    for central_ring in benzene_rings:
        attachments = []  # list of tuples: (central_atom_position_index, peripheral_ring_id)
        # For each atom (with its index in the ring order) in the candidate central ring:
        for pos, atom_idx in enumerate(central_ring):
            atom = mol.GetAtomWithIdx(atom_idx)
            # Look at each neighbor that is NOT in the central ring.
            for nb in atom.GetNeighbors():
                nb_idx = nb.GetIdx()
                if nb_idx in central_ring:
                    continue
                # We pass the set of indices of the central ring as "excluded" (converted to a set)
                peripheral_ring_id = dfs_find_peripheral(nb_idx, excluded=set(central_ring))
                if peripheral_ring_id is not None:
                    attachments.append((pos, peripheral_ring_id))
        # For a valid para-terphenyl, we want at least two distinct peripheral benzene rings attached.
        if len(attachments) < 2:
            continue

        # Check pairs of attachments for (a) distinct peripheral rings and (b) para connection on the central ring.
        for i in range(len(attachments)):
            for j in range(i+1, len(attachments)):
                pos_i, per_id_i = attachments[i]
                pos_j, per_id_j = attachments[j]
                # They must be attached to two distinct benzene rings.
                if per_id_i == per_id_j:
                    continue
                # Compute the circular distance along the six-membered ring.
                diff = abs(pos_i - pos_j)
                circ_dist = min(diff, 6 - diff)
                if circ_dist == 3:
                    return True, ("Para-terphenyl core found: central benzene ring substituted at para positions " +
                                  "with two distinct benzene rings (allowing for a short linker up to 2 bonds).")
    
    return False, "Para-terphenyl core not found in the molecule."

# Debug/test section (can be removed or commented out in production)
if __name__ == "__main__":
    # Some illustrative test cases:
    test_cases = {
        "1,4-diphenylbenzene": "c1ccc(cc1)-c2ccc(cc2)-c3ccccc3",
        "Benzene": "c1ccccc1",
        # The provided examples (only one displayed here; more can be added for testing):
        "3-Methoxyterprenin": "O(C1=C(OC)C=C(C2=C(O)C(OC)=C(C3=CC=C(O)C=C3)C=C2OC)C=C1)CC=C(C)C",
        "Kynapcin-12": "O=C(OC1=C(O)C(=C(OC(=O)C)C(=C1C2=CC=C(O)C=C2)O)C3=CC=C(O)C=C3)C"
    }
    for name, smi in test_cases.items():
        result, reason = is_para_terphenyl(smi)
        print(f"Test: {name}\n SMILES: {smi}\n Result: {result} | Reason: {reason}\n")