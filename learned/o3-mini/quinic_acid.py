"""
Classifies: CHEBI:26493 quinic acid
"""
"""
Classifies: Quinic acid – a cyclitol carboxylic acid derivative.
A quinic acid derivative has a saturated cyclohexane ring (all carbons,
with only single bonds) that displays exactly one exocyclic carboxyl (–C(=O)O–)
and several oxygen substituents (typically four, whether free hydroxyls or esterified).
This is a heuristic classifier.
"""
from rdkit import Chem

def is_quinic_acid(smiles: str):
    """
    Determines if a molecule is a quinic acid (cyclitol carboxylic acid derivative)
    based on its SMILES string. The algorithm looks for a saturated six-membered carbon
    ring (cyclohexane) that has exactly one exocyclic carboxyl group (-C(=O)O-) and at least
    four oxygen substituents attached via single bonds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is recognized as a quinic acid derivative, False otherwise.
        str: Reason describing the outcome.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    def is_carboxyl_group(carbon):
        """
        Checks if the given carbon atom is part of a carboxyl (or ester carboxyl) group.
        A proper carboxyl carbon should be sp2 (often it will make one double bond to oxygen
        and one single bond to oxygen) regardless of additional substituents.
        """
        if carbon.GetAtomicNum() != 6:
            return False
        oxy_double = 0
        oxy_single = 0
        # Loop over neighbors to check bonded oxygen count and bond types.
        for nbr in carbon.GetNeighbors():
            if nbr.GetAtomicNum() == 8:
                bond = mol.GetBondBetweenAtoms(carbon.GetIdx(), nbr.GetIdx())
                if bond is None:
                    continue
                if bond.GetBondType() == Chem.BondType.DOUBLE:
                    oxy_double += 1
                elif bond.GetBondType() == Chem.BondType.SINGLE:
                    oxy_single += 1
        # Heuristic: require at least one double- and one single-bonded oxygen.
        if oxy_double >= 1 and oxy_single >= 1:
            return True
        return False

    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()

    # Examine each ring for the desired cyclohexane substructure.
    for ring in atom_rings:
        if len(ring) != 6:
            continue  # Only consider six-membered rings.
        # Ensure all atoms in the ring are carbons.
        if not all(mol.GetAtomWithIdx(idx).GetAtomicNum() == 6 for idx in ring):
            continue
        # Ensure the ring bonds are all single (i.e. fully saturated, non-aromatic)
        ring_is_saturated = True
        n_ring = len(ring)
        for i in range(n_ring):
            a_idx = ring[i]
            b_idx = ring[(i+1) % n_ring]  # Next atom, wrap around.
            bond = mol.GetBondBetweenAtoms(a_idx, b_idx)
            if bond is None:
                ring_is_saturated = False
                break
            # Also check that the bond is not aromatic
            if bond.GetIsAromatic() or bond.GetBondType() != Chem.BondType.SINGLE:
                ring_is_saturated = False
                break
        if not ring_is_saturated:
            continue

        # Now count exocyclic substituents:
        carboxyl_count = 0
        oxy_substituent_count = 0
        # Look over each atom of the ring and its external neighbors.
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            for nbr in atom.GetNeighbors():
                # Only consider substituents not in the ring.
                if nbr.GetIdx() in ring:
                    continue
                bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
                if nbr.GetAtomicNum() == 8 and bond.GetBondType() == Chem.BondType.SINGLE:
                    # Count as an oxygen substituent (e.g. as in an –OH or part of an ester linkage)
                    oxy_substituent_count += 1
                elif nbr.GetAtomicNum() == 6:
                    # Check if the neighbor carbon qualifies as a carboxyl group.
                    if is_carboxyl_group(nbr):
                        carboxyl_count += 1
        # Heuristic: for quinic acid, require exactly one carboxyl and at least four oxygen substituents.
        if carboxyl_count == 1 and oxy_substituent_count >= 4:
            return True, (f"Contains a saturated cyclohexane ring with {carboxyl_count} carboxyl group "
                          f"and {oxy_substituent_count} oxygen substituents, consistent with a cyclitol carboxylic acid derivative.")
    return False, "No qualifying saturated cyclohexane ring with one carboxyl group and at least four oxygen substituents found"

# (Optional testing code)
# smiles_examples = [
#     "O[C@H]1C[C@@](O)(C[C@H](O)[C@H]1O)C(O)=O",  # (+)-quinic acid
#     "COc1cc(ccc1O)\\C=C\\C(=O)O[C@H]1[C@H](O)C[C@@](O)(C[C@H]1O)C(O)=O",  # 4-O-feruloyl-D-quinic acid
# ]
# for smi in smiles_examples:
#     result, reason = is_quinic_acid(smi)
#     print(smi, "->", result, reason)