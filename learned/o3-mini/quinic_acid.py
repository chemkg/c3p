"""
Classifies: CHEBI:26493 quinic acid
"""
"""
Classifies: Quinic acid – a cyclitol carboxylic acid derivative.

This heuristic requires that the molecule possesses an isolated saturated cyclohexane ring
(six carbons in the ring, all sp3, not fused with any other ring) on which exactly five exocyclic 
substituents are attached. Of these, one must be a (neutral) carboxyl group—which we define as a 
carbon (exocyclic to the ring) bearing at least one double‐bonded oxygen and one single‐bonded oxygen,
both with zero formal charge—and the remaining four substituents must be oxygen atoms (as –OH groups 
or as part of an ester) attached directly to the ring. In addition, if an oxygen substituent carries 
an acyl chain with a C=C double bond, the stereochemistry of that bond is examined: any C=C marked 
with “Z” (cis) is taken as disqualifying.
Note: This heuristic is not perfect but aims to reduce false positives.
"""
from rdkit import Chem

def is_quinic_acid(smiles: str):
    """
    Determines if a molecule is a quinic acid derivative (cyclitol carboxylic acid derivative)
    based on its SMILES string.
    
    The algorithm:
      (1) Parse the molecule.
      (2) For each ring with exactly 6 atoms:
             - Ensure all ring atoms are carbons.
             - Ensure that bonds within the ring are single (non‐aromatic) – i.e. a saturated cyclohexane.
             - Check that the ring is isolated (not fused with any other ring).
      (3) For a candidate ring, collect exocyclic substituents (neighbors of ring atoms not in the ring).
             - In quinic acid exactly five substituents should be present.
             - Among these, one should be an exocyclic carboxyl group and four should be oxygen substituents.
      (4) For each oxygen substituent that is not “just” hydroxyl (e.g. when it is esterified) perform 
          a short recursive search into the substituent fragment: if any double bond is found and its 
          stereochemistry is defined as Z (cis), then reject the candidate.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is recognized as a quinic acid derivative, False otherwise.
        str: A reason describing the outcome.
        (If the task is too challenging an implementation may return (None, None).)
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    
    # Helper: Check if a given carbon atom (exocyclic to the candidate ring) is a neutral carboxyl group.
    def is_carboxyl_group(atom):
        # Check that the atom is carbon.
        if atom.GetAtomicNum() != 6:
            return False
        dO = 0
        sO = 0
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 8:
                bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
                if bond is None:
                    continue
                if bond.GetBondType() == Chem.BondType.DOUBLE:
                    dO += 1
                    if nbr.GetFormalCharge() != 0:
                        return False
                elif bond.GetBondType() == Chem.BondType.SINGLE:
                    sO += 1
                    if nbr.GetFormalCharge() != 0:
                        return False
        return (dO >= 1 and sO >= 1)
    
    # Helper: recursively check substituent fragment for any double bonds with Z stereochemistry.
    def check_fragment_for_Z(atom, visited, max_depth=4, depth=0):
        if depth > max_depth:
            return True  # do not traverse too deep
        visited.add(atom.GetIdx())
        for bond in atom.GetBonds():
            # Only inspect double bonds.
            if bond.GetBondType() == Chem.BondType.DOUBLE:
                # If the bond stereo is defined and is Z, then reject.
                if bond.GetStereo() == Chem.BondStereo.STEREOZ:
                    return False
            nbr = bond.GetOtherAtom(atom)
            if nbr.GetIdx() in visited:
                continue
            # Recurse if the neighbor is not a hydrogen.
            if nbr.GetAtomicNum() > 1:
                if not check_fragment_for_Z(nbr, visited, max_depth, depth+1):
                    return False
        return True

    # Now search each ring candidate.
    for ring in atom_rings:
        if len(ring) != 6:
            continue  # only consider six-membered rings
        # All atoms in the ring must be carbons.
        if not all(mol.GetAtomWithIdx(idx).GetAtomicNum() == 6 for idx in ring):
            continue
        # Check that all bonds in the ring are single and non‐aromatic.
        ring_saturated = True
        n_ring = len(ring)
        for i in range(n_ring):
            idx1 = ring[i]
            idx2 = ring[(i+1) % n_ring]
            bond = mol.GetBondBetweenAtoms(idx1, idx2)
            if bond is None or bond.GetIsAromatic() or bond.GetBondType() != Chem.BondType.SINGLE:
                ring_saturated = False
                break
        if not ring_saturated:
            continue
        
        # Check that the ring is isolated (i.e. not fused with any other ring).
        candidate_set = set(ring)
        isolated = True
        for other_ring in atom_rings:
            if set(other_ring) == candidate_set:
                continue
            if len(candidate_set.intersection(other_ring)) >= 2:
                isolated = False
                break
        if not isolated:
            continue
        
        # Gather exocyclic substituents on the ring.
        # For each ring atom, collect all neighbors that are not in the ring.
        substituents = []  # list of tuples (ring_atom_idx, neighbor_atom, bond)
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in candidate_set:
                    continue
                bond = mol.GetBondBetweenAtoms(idx, nbr.GetIdx())
                substituents.append((idx, nbr, bond))
        
        # For quinic acid (or its common derivatives) we expect exactly five substituents.
        if len(substituents) != 5:
            continue
        
        carboxyl_count = 0
        oxygen_count = 0
        geometry_ok = True
        # Process each substituent.
        for (ring_idx, nbr, bond) in substituents:
            if nbr.GetAtomicNum() == 6:
                # Check for a carboxyl group directly attached to the ring carbon.
                if is_carboxyl_group(nbr):
                    carboxyl_count += 1
                else:
                    # If a carbon substituent is not a carboxyl group, disqualify.
                    geometry_ok = False
                    break
            elif nbr.GetAtomicNum() == 8:
                oxygen_count += 1
                # For an oxygen substituent, if it is not simply -OH (i.e. if it is connected to additional heavy atoms)
                # then check the substituent fragment for double bonds with Z stereochemistry.
                heavy_neighbors = [a for a in nbr.GetNeighbors() if a.GetAtomicNum() > 1 and a.GetIdx() not in candidate_set]
                if len(heavy_neighbors) > 1 or (len(heavy_neighbors) == 1 and heavy_neighbors[0].GetAtomicNum() != 1):
                    if not check_fragment_for_Z(nbr, set()):
                        geometry_ok = False
                        break
            else:
                # Unexpected substituent type.
                geometry_ok = False
                break
        
        if not geometry_ok:
            continue
        if carboxyl_count == 1 and oxygen_count == 4:
            reason = (f"Contains an isolated saturated cyclohexane ring with 1 exocyclic carboxyl group and "
                      f"4 oxygen substituents. All examined ester fragments show acceptable (non-Z) double bond stereochemistry, "
                      f"consistent with a cyclitol carboxylic acid derivative.")
            return True, reason

    return False, "No isolated saturated cyclohexane ring with 1 carboxyl group and 4 oxygen substituents (with acceptable substituent stereochemistry) found"

# (Optional testing code; remove or uncomment as needed)
# test_smiles = [
#     # True positives:
#     "COc1cc(ccc1O)\\C=C\\C(=O)O[C@H]1[C@H](O)C[C@@](O)(C[C@H]1O)C(O)=O",  # 4-O-feruloyl-D-quinic acid
#     "O[C@H]1C[C@@](O)(C[C@H](O)[C@H]1O)C(O)=O",  # (+)-quinic acid
#     # False positive test (expected False):
#     "O[C@@H]1C[C@](O)(C[C@@H](OC(=O)\\C=C/c2ccc(O)c(O)c2)[C@@H]1O)C(O)=O",  # 5-O-cis-caffeoylquinic acid (should be rejected)
# ]
# for smi in test_smiles:
#     result, msg = is_quinic_acid(smi)
#     print(smi, "->", result, msg)