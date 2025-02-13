"""
Classifies: CHEBI:26493 quinic acid
"""
"""
Classifies: Quinic acid – a cyclitol carboxylic acid derivative.
A quinic acid derivative is defined as having an isolated (i.e. non-fused)
saturated cyclohexane ring (all carbons, all single bonds) that displays exactly
one exocyclic carboxyl group (with neutral oxygen atoms) and at least four oxygen substituents.
This heuristic attempts to avoid false positives from fused-ring molecules.
"""
from rdkit import Chem

def is_quinic_acid(smiles: str):
    """
    Determines if a molecule is a quinic acid derivative (cyclitol carboxylic acid derivative)
    based on its SMILES string.
    The algorithm:
      (1) Searches for a saturated cyclohexane (six-carbon) ring.
      (2) Ensures that the ring is not fused with any other ring (i.e. it is isolated).
      (3) Looks for exactly one exocyclic carboxyl group (as defined by a carbon
          that is sp2 and bound to one double-bonded oxygen and one single-bonded oxygen, both neutral).
      (4) Requires at least four oxygen substituents from the ring carbons.
    Args:
        smiles (str): SMILES string of the molecule.
    Returns:
        bool: True if the molecule is recognized as a quinic acid derivative, False otherwise.
        str: A reason describing the outcome.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    
    # Helper function to check if a carbon is part of a neutral carboxyl group.
    def is_carboxyl_group(carbon):
        """
        Confirms whether a carbon atom is part of a carboxyl group.
        It must be a carbon (atomic number 6) with at least one double-bonded oxygen and one single-bonded oxygen,
        and the involved oxygens must have a formal charge of zero.
        """
        if carbon.GetAtomicNum() != 6:
            return False
        oxy_double = []
        oxy_single = []
        for nbr in carbon.GetNeighbors():
            if nbr.GetAtomicNum() == 8:
                bond = mol.GetBondBetweenAtoms(carbon.GetIdx(), nbr.GetIdx())
                if bond is None:
                    continue
                if bond.GetBondType() == Chem.BondType.DOUBLE:
                    oxy_double.append(nbr)
                elif bond.GetBondType() == Chem.BondType.SINGLE:
                    oxy_single.append(nbr)
        # Require at least one double and one single bond to oxygen and check formal charges are zero.
        if oxy_double and oxy_single:
            if all(o.GetFormalCharge() == 0 for o in oxy_double + oxy_single):
                return True
        return False

    # Check each ring to find an isolated cyclohexane ring.
    for ring in atom_rings:
        if len(ring) != 6:
            continue  # Only consider six-membered rings.
        
        # Ensure all atoms in the ring are carbons.
        if not all(mol.GetAtomWithIdx(idx).GetAtomicNum() == 6 for idx in ring):
            continue
            
        # Ensure that all bonds within the ring are single and not aromatic.
        ring_is_saturated = True
        n_ring = len(ring)
        for i in range(n_ring):
            idx1 = ring[i]
            idx2 = ring[(i+1) % n_ring]
            bond = mol.GetBondBetweenAtoms(idx1, idx2)
            if bond is None or bond.GetIsAromatic() or bond.GetBondType() != Chem.BondType.SINGLE:
                ring_is_saturated = False
                break
        if not ring_is_saturated:
            continue
            
        # Now check that the cyclohexane ring is isolated (not fused with another ring).
        # For every other ring in the molecule, if it shares 2 or more atoms with our candidate ring, then discard.
        isolated = True
        candidate_set = set(ring)
        for other_ring in atom_rings:
            if set(other_ring) == candidate_set:
                continue  # this is our candidate ring itself
            # If the intersection is 2 or more atoms, assume the candidate ring is fused.
            if len(candidate_set.intersection(other_ring)) >= 2:
                isolated = False
                break
        if not isolated:
            continue

        # Now count substituents on the candidate ring.
        carboxyl_count = 0
        oxy_substituent_count = 0
        # For each carbon in the ring, examine neighbors not in the ring.
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in candidate_set:
                    continue
                bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
                # If the neighbor is oxygen and connected by a single bond, count it.
                if nbr.GetAtomicNum() == 8 and bond.GetBondType() == Chem.BondType.SINGLE:
                    oxy_substituent_count += 1
                # If the neighbor is carbon, test if it is a carboxyl group.
                elif nbr.GetAtomicNum() == 6:
                    if is_carboxyl_group(nbr):
                        # Also ensure that the attached carboxyl group is neutral.
                        carboxyl_count += 1
        # For a quinic acid derivative, we require exactly one carboxyl group and at least four oxygens.
        if carboxyl_count == 1 and oxy_substituent_count >= 4:
            reason = (f"Contains an isolated saturated cyclohexane ring with {carboxyl_count} carboxyl group and "
                      f"{oxy_substituent_count} oxygen substituents, consistent with a cyclitol carboxylic acid derivative.")
            return True, reason
    
    return False, "No isolated saturated cyclohexane ring with one (neutral) carboxyl group and at least four oxygen substituents found"

# (Optional testing code; remove or uncomment as needed)
# examples = [
#     "COc1cc(ccc1O)\\C=C\\C(=O)O[C@H]1[C@H](O)C[C@@](O)(C[C@H]1O)C(O)=O",  # 4-O-feruloyl-D-quinic acid, correct
#     "O[C@@H]1C[C@@](O)(C[C@@H](O)[C@H]1O)C(O)=O",  # (+)-quinic acid, correct
#     "O[C@@H]1C[C@@](O)(C[C@@H](O)[C@H]1O)C([O-])=O",  # (-)-quinate, should be rejected
# ]
# for smi in examples:
#     result, msg = is_quinic_acid(smi)
#     print(smi, "->", result, msg)