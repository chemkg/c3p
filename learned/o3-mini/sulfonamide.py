"""
Classifies: CHEBI:35358 sulfonamide
"""
"""
Classifies: Sulfonamide
Definition: An amide of a sulfonic acid, RS(=O)(=O)NR'2.
This module checks if a molecule contains a valid sulfonamide group.
The approach is to first flag candidate S(=O)(=O)-N motifs using a SMARTS pattern,
then for each candidate verify that:
  1. The sulfur atom has exactly two double bonds to oxygen.
  2. It has exactly one bond to a nitrogen.
  3. It has exactly one additional substituent (typically R) so that its total degree is 4.
If a candidate passes those checks we classify the molecule as a sulfonamide.
"""

from rdkit import Chem

def is_sulfonamide(smiles: str):
    """
    Determines if a molecule is a sulfonamide based on its SMILES string.
    A sulfonamide is defined as RS(=O)(=O)NR'2.
    
    The algorithm first uses a SMARTS pattern to pick up any potential S(=O)(=O)-N
    connectivity. Then for each candidate it inspects the sulfur atom’s bonds to
    ensure that it has exactly two double bonds to oxygen, one bond to a nitrogen,
    and a fourth bond to an R group.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a sulfonamide, False otherwise.
        str: Reason for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # The basic SMARTS pattern to catch an S(=O)(=O)-N fragment.
    sulfonamide_smarts = "[S](=O)(=O)-[N]"
    pattern = Chem.MolFromSmarts(sulfonamide_smarts)
    matches = mol.GetSubstructMatches(pattern)
    if not matches:
        return False, "No sulfonamide group (S(=O)(=O)-N) found in the molecule"

    # Now, check each match in detail.
    for match in matches:
        # By design the SMARTS returns two indices: first is the sulfur, second is the nitrogen.
        s_idx, n_idx = match[0], match[1]
        s_atom = mol.GetAtomWithIdx(s_idx)
        n_atom = mol.GetAtomWithIdx(n_idx)

        # Check the number of bonds from sulfur that are double bonds to oxygen.
        dO_count = 0
        for bond in s_atom.GetBonds():
            neighbor = bond.GetOtherAtom(s_atom)
            if neighbor.GetAtomicNum() == 8 and bond.GetBondTypeAsDouble() == 2:
                dO_count += 1
        if dO_count != 2:
            # This candidate S does not have exactly two S=O bonds.
            continue

        # Count how many bonds from sulfur go to nitrogen.
        n_bond_count = sum(1 for bond in s_atom.GetBonds() if bond.GetOtherAtom(s_atom).GetAtomicNum() == 7)
        if n_bond_count != 1:
            # There should be exactly one bond from S to nitrogen.
            continue

        # Ensure that the sulfur has exactly one more (non-O, non-N) substituent.
        # In a typical sulfonamide, S is tetravalent: two oxygens (double bonds), one N, and one R group.
        # Note: GetDegree() counts each bonded neighbor once.
        if s_atom.GetDegree() != 4:
            continue

        # Identify the substituents of S that are NOT the nitrogen and not oxygens.
        other_substituents = [nbr for nbr in s_atom.GetNeighbors() if nbr.GetAtomicNum() not in [7, 8]]
        if len(other_substituents) != 1:
            continue

        # If all conditions are met, we can accept this sulfonamide candidate:
        return True, "Molecule contains a sulfonamide group: RS(=O)(=O)NR'2"

    # If no candidate passed all filters, then we return failure.
    return False, "No valid sulfonamide group found in the molecule"

# Uncomment below to run some quick tests.
# test_smiles = [
#     "Brc1ccc(c2ccccc12)S(=O)(=O)NCc1ccccn1",  # pyrabactin, expected True
#     "CC1=CC=C(C=C1)S(=O)(=O)N2CC3(C2)CN([C@@H](C4=C3C5=C(N4)C=C(C=C5)OC)CO)CC6=CC=CC=C6",  # additional example
# ]
# for sm in test_smiles:
#     result, reason = is_sulfonamide(sm)
#     print(f"SMILES: {sm}\nResult: {result}\nReason: {reason}\n")