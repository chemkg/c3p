"""
Classifies: CHEBI:50753 isoflavonoid
"""
"""
Classifies: isoflavonoid
Definition: Any 1-benzopyran with an aryl substituent at position 3.
The algorithm first locates fused rings that match a benzene (6 aromatic carbons) 
and a pyran (6 atoms: exactly one oxygen, five carbons). Then it identifies the atoms 
of the pyran that are not shared with the benzene (the “free” pyran atoms) and 
further restricts to only those free atoms that are adjacent to at least one fused atom.
Finally, among these candidates the algorithm requires that exactly one has a single‐bond 
neighbor that is (a) aromatic and (b) belongs to a benzene ring (6 atoms, all aromatic carbons). 
If so, the molecule is classified as an isoflavonoid.
Note: This heuristic may not catch every isoflavonoid but is designed to improve
on the previously observed false‐positives and false‐negatives.
"""

from rdkit import Chem

def is_isoflavonoid(smiles: str):
    """
    Determines if a molecule is an isoflavonoid based on its SMILES string.
    An isoflavonoid (for our purpose) is defined as a 1-benzopyran core 
    (a benzene ring fused to a pyran ring, where the pyran ring has one oxygen and five carbons)
    with an aryl (benzene) substituent attached at the free (non‐fused) carbon of the pyran 
    (expected to be position 3).
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as an isoflavonoid, False otherwise.
        str: Explanation of the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    ring_info = mol.GetRingInfo().AtomRings()
    if not ring_info:
        return False, "No rings found in molecule"
    
    # Helper: check if a ring is benzene: 6 atoms, all carbons and aromatic.
    def is_benzene(ring):
        if len(ring) != 6:
            return False
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() != 6 or not atom.GetIsAromatic():
                return False
        return True

    # Helper: check if a ring is pyran: 6 atoms with exactly one oxygen and five carbons.
    def is_pyran(ring):
        if len(ring) != 6:
            return False
        o_count = 0
        c_count = 0
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() == 8:
                o_count += 1
            elif atom.GetAtomicNum() == 6:
                c_count += 1
            else:
                return False
        return (o_count == 1 and c_count == 5)
    
    # We want to try each pair of rings that share exactly 2 atoms.
    # One must be benzene and the other pyran.
    candidate_found = False
    reason = "No valid benzopyran core with a 3-aryl substituent was detected"
    for i in range(len(ring_info)):
        for j in range(i+1, len(ring_info)):
            ring1 = ring_info[i]
            ring2 = ring_info[j]
            shared = set(ring1).intersection(set(ring2))
            if len(shared) != 2:
                continue  # not fused by exactly two atoms
            # Determine which is benzene and which is pyran:
            benzene_ring = None
            pyran_ring = None
            if is_benzene(ring1) and is_pyran(ring2):
                benzene_ring = set(ring1)
                pyran_ring = set(ring2)
            elif is_benzene(ring2) and is_pyran(ring1):
                benzene_ring = set(ring2)
                pyran_ring = set(ring1)
            else:
                continue
            
            # Get free (non-fused) atoms of the pyran ring.
            free_pyran = list(pyran_ring - shared)
            # To restrict the candidate to the expected (position 3) substitution,
            # we require that the candidate free pyran atom is (a) a carbon,
            # (b) is adjacent in the pyran ring to at least one fused (shared) atom,
            # and (c) carries an aryl substituent (defined below).
            candidate_substitutions = []
            for idx in free_pyran:
                atom = mol.GetAtomWithIdx(idx)
                if atom.GetAtomicNum() != 6:
                    continue  # substitution expected on carbon only
                # Check if this free atom is adjacent (via pyran bond) to a fused atom.
                adjacent_fused = False
                for nbr in atom.GetNeighbors():
                    if nbr.GetIdx() in shared:
                        adjacent_fused = True
                        break
                if not adjacent_fused:
                    continue
                # Now look for neighbors outside the pyran ring that could be the aryl substituent.
                for nbr in atom.GetNeighbors():
                    if nbr.GetIdx() in pyran_ring:
                        continue  # ignore atoms in core
                    # Only consider single bonds (as substituents are usually attached by a single bond)
                    bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
                    if bond is None or bond.GetBondTypeAsDouble() != 1.0:
                        continue
                    # Check that neighbor is aromatic.
                    if not nbr.GetIsAromatic():
                        continue
                    # Now confirm that the neighbor is part of a benzene ring.
                    found_benzene = False
                    for ring in ring_info:
                        if len(ring) == 6 and nbr.GetIdx() in ring:
                            # All atoms in this ring must be aromatic carbons.
                            if all(mol.GetAtomWithIdx(a).GetAtomicNum() == 6 and mol.GetAtomWithIdx(a).GetIsAromatic() for a in ring):
                                found_benzene = True
                                break
                    if found_benzene:
                        candidate_substitutions.append((atom.GetIdx(), nbr.GetIdx()))
                        # We break after one validated neighbor per candidate atom.
                        break

            # We now require that exactly one free pyran atom was found to have a valid aryl substituent.
            if len(candidate_substitutions) == 1:
                return (True, 
                        ("Found benzopyran core (benzene fused with a pyran ring having 1 O and 5 C) with an aryl substituent "
                         "attached at a free pyran carbon (likely position 3)."))
            elif len(candidate_substitutions) > 1:
                return (False, 
                        f"Multiple aryl substituents ({len(candidate_substitutions)}) found on the pyran core; expected exactly one.")
            # Otherwise, continue checking other fused ring pairs.
    
    return (False, reason)

# Example usage (for testing the function):
if __name__ == "__main__":
    test_smiles = [
        "Oc1ccc(cc1)-c1coc2cc(O)ccc2c1=O",  # daidzein, isoflavonoid; expected True
        "CC(=O)O",                        # acetic acid, negative control; expected False
        "O=c1c2c(O)cc(O)cc2oc2cc(O)ccc12",  # example that may be borderline
    ]
    for smi in test_smiles:
        classification, reason = is_isoflavonoid(smi)
        print(f"SMILES: {smi}\nClassified as isoflavonoid: {classification}\nReason: {reason}\n")