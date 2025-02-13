"""
Classifies: CHEBI:35358 sulfonamide
"""
"""
Classifies: Sulfonamide
Definition: An amide of a sulfonic acid, RS(=O)2NR'2.
This module checks if a molecule contains a valid sulfonamide group.
The approach is to first flag a candidate S(=O)(=O)-N motif via SMARTS,
then verify that the sulfur atom has exactly two double-bonded oxygens,
and an overall connectivity consistent with a sulfonamide:
  • The sulfur normally has four bonds: one to R, two double bonds to O, 
    and one single bond to N.
This should help reduce false positives (where the moiety is embedded elsewhere)
and false negatives.
"""

from rdkit import Chem

def is_sulfonamide(smiles: str):
    """
    Determines if a molecule is a sulfonamide based on its SMILES string.
    A sulfonamide is defined as RS(=O)(=O)NR'2.
    
    The algorithm first finds a candidate group using a SMARTS pattern that matches a 
    sulfur atom double-bonded to two oxygens and singly bonded to a nitrogen.
    It then checks that the sulfur atom has exactly two double bonds to oxygen and 
    an overall degree consistent with an RS(=O)(=O)NR'2 moiety.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a sulfonamide, False otherwise.
        str: Reason for the classification decision.
    """
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern to identify an S(=O)(=O)-N connection.
    # (This flag is intentionally less restrictive so that subsequent checks can filter.)
    sulfonamide_smarts = "[S](=O)(=O)-[N]"
    pattern = Chem.MolFromSmarts(sulfonamide_smarts)
    matches = mol.GetSubstructMatches(pattern)
    if not matches:
        return False, "No sulfonamide group (S(=O)(=O)-N) found in the molecule"
    
    # For every match, verify the connectivity and bond orders around the sulfur atom.
    # We assume the pattern order: match[0] is the sulfur, match[1] is the nitrogen.
    for match in matches:
        s_idx = match[0]
        n_idx = match[1]
        s_atom = mol.GetAtomWithIdx(s_idx)
        # Count bonds from sulfur that are double bonds to oxygen.
        dO_count = 0
        for bond in s_atom.GetBonds():
            # Check if neighbor is oxygen with a double bond.
            neighbor = bond.GetOtherAtom(s_atom)
            # Use bond.GetBondTypeAsDouble() which returns 2.0 if double.
            if neighbor.GetAtomicNum() == 8 and bond.GetBondTypeAsDouble() == 2:
                dO_count += 1
        
        # In a typical sulfonamide RS(=O)(=O)NR'2, the sulfur should be tetravalent
        # with exactly two double-bonded oxygens, one bond to nitrogen and one to an R substituent.
        if dO_count != 2:
            continue  # This match does not meet our S=O criteria.
        
        # To be extra sure, we check that sulfur's total degree is 4.
        # (S is bonded to two O (double bonds count as one connection each),
        # one N from the sulfonamide and one carbon or other substituent for R.)
        if s_atom.GetDegree() != 4:
            continue
        
        # At this stage the candidate S has two double-bonded oxygens.
        # We also verify that the bonded nitrogen (n_idx) exists.
        n_atom = mol.GetAtomWithIdx(n_idx)
        # Do not enforce a rigid count on N neighbors here because the nitrogen in many
        # sulfonamides can have additional bonds (e.g., aromatic rings or extra substituents).
        # If we reached here, then the pattern and connectivity look like a sulfonamide.
        return True, "Molecule contains a sulfonamide group: RS(=O)(=O)NR'2"
    
    # If no valid match is found after checking all candidates, then we fail.
    return False, "No valid sulfonamide group found in the molecule"

# Uncomment below to test a few SMILES examples.
# test_smiles = [
#     "Brc1ccc(c2ccccc12)S(=O)(=O)NCc1ccccn1",  # pyrabactin, expected True
#     "CC1=CC=C(C=C1)S(=O)(=O)N2CC3(C2)CN([C@@H](C4=C3C5=C(N4)C=C(C=C5)OC)CO)CC6=CC=CC=C6",  # another example
# ]
# for sm in test_smiles:
#     result, reason = is_sulfonamide(sm)
#     print(f"SMILES: {sm}\nResult: {result}\nReason: {reason}\n")