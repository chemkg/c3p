"""
Classifies: CHEBI:26493 quinic acid
"""
"""
Classifies: Quinic acid – a cyclitol carboxylic acid.
A quinic acid derivative has a cyclohexane ring (non‐aromatic, all C) that displays
an exocyclic carboxyl (–C(=O)O–) and several oxygen substituents (which in the parent are free –OH).
Note that many quinates exist where the –OH groups are esterified.
This is just one heuristic classifier.
"""
from rdkit import Chem

def is_quinic_acid(smiles: str):
    """
    Determines if a molecule is a quinic acid (cyclitol carboxylic acid) derivative based on its SMILES string.
    The algorithm looks for a six‐membered non‐aromatic (cyclohexane) ring made solely of carbon atoms,
    in which at least one ring carbon is substituted with a carboxyl group –C(=O)O– and the ring exhibits
    several oxygen substituents (as free hydroxyls or ester groups).

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

    # Helper function: Check if an atom (assumed to be carbon) is a carboxyl (acid/ester) carbon.
    def is_carboxyl_group(carbon):
        # Carboxyl (or its ester derivative) has a carbon bound to two oxygens:
        # one by a double bond and one by a single bond.
        if carbon.GetAtomicNum() != 6:
            return False
        # We require exactly two oxygen neighbors (aside from the connection to the ring)
        oxy_double = 0
        oxy_single = 0
        # Loop over all neighbors of the carbon candidate.
        for nbr in carbon.GetNeighbors():
            if nbr.GetAtomicNum() == 8:
                bond = mol.GetBondBetweenAtoms(carbon.GetIdx(), nbr.GetIdx())
                if bond is None:
                    continue
                if bond.GetBondType() == Chem.BondType.DOUBLE:
                    oxy_double += 1
                elif bond.GetBondType() == Chem.BondType.SINGLE:
                    oxy_single += 1
        # A proper carboxyl (or ester carboxyl) carbon should have at least one double-bonded oxygen
        # and one single-bonded oxygen.
        if oxy_double >= 1 and oxy_single >= 1:
            return True
        return False

    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()

    # Examine each ring for the desired substructure: cyclohexane
    for ring in atom_rings:
        if len(ring) != 6:
            continue  # Only consider 6-membered rings
        # Check that all atoms in this ring are carbons (i.e. cyclohexane)
        if not all(mol.GetAtomWithIdx(idx).GetAtomicNum() == 6 for idx in ring):
            continue

        carboxyl_found = False
        oxygen_substituent_count = 0
        # For each atom in the ring, look at neighbors not in the ring
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring:
                    continue  # ignore atoms in the ring
                bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
                # If the neighbor is an oxygen and the bond is a single bond,
                # we count it as a hydroxyl (or an oxygen substituent)
                if nbr.GetAtomicNum() == 8 and bond.GetBondType() == Chem.BondType.SINGLE:
                    oxygen_substituent_count += 1
                # If the neighbor is a carbon, check if it is a carboxyl group.
                elif nbr.GetAtomicNum() == 6:
                    # Check if this carbon acts as a carboxyl group
                    if is_carboxyl_group(nbr):
                        carboxyl_found = True
        # Heuristic: in quinic acid the cyclitol ring has one carboxyl group and typically 3 or more oxygen substituents.
        if carboxyl_found and oxygen_substituent_count >= 3:
            return True, ("Contains a cyclohexane ring with a carboxyl substituent and "
                          f"{oxygen_substituent_count} oxygen substituents, consistent with a cyclitol carboxylic acid.")

    return False, "No qualifying cyclohexane ring with the required carboxyl and oxygen substituents found"

# Testing (optional): Uncomment below lines to test on one example
# smiles_examples = [
#     "O[C@H]1C[C@@](O)(C[C@H](O)[C@H]1O)C(O)=O",  # (+)-quinic acid
#     "COc1cc(ccc1O)\\C=C\\C(=O)O[C@H]1[C@H](O)C[C@@](O)(C[C@H]1O)C(O)=O",  # 4-O-feruloyl-D-quinic acid
# ]
# for smi in smiles_examples:
#     result, reason = is_quinic_acid(smi)
#     print(smi, "->", result, reason)