"""
Classifies: CHEBI:22307 aldoxime
"""
"""
Classifies: Aldoxime (Oximes of aldehydes RCH=NOH)
An aldoxime is formed from an aldehyde (R–CHO) by conversion of the carbonyl group to an oxime group (R–CH=N–OH).
Thus, the structure must contain an aldoxime functional group with:
  - A carbon atom with exactly one hydrogen (the aldehyde carbon) involved in a double bond,
  - A nitrogen atom with exactly one hydrogen,
  - And an oxygen atom bearing one hydrogen (–OH) attached via a single bond to the nitrogen.
This script uses a SMARTS pattern to define that substructure.
"""

from rdkit import Chem

def is_aldoxime(smiles: str):
    """
    Determines if a molecule is an aldoxime (RCH=NOH) based on its SMILES string.
    This function converts the SMILES into an RDKit molecule (with explicit hydrogens added)
    and then searches for the aldoxime functional group using a SMARTS pattern.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        (bool, str): A tuple where the first element is True if the aldoxime group is found,
                     and False otherwise; the second element is a string providing a reason.
    """
    # Convert the SMILES string to a molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens so that hydrogen counts are correctly represented.
    mol = Chem.AddHs(mol)
    
    # Define a SMARTS pattern for an aldoxime group: R–CH=N–OH.
    # The pattern [C;H1] ensures the carbon has exactly one hydrogen,
    # [N;H1] ensures the nitrogen carries one hydrogen,
    # and [O;H1] ensures the oxygen is present as an -OH.
    aldoxime_smarts = "[C;H1]=[N;H1][O;H1]"
    pattern = Chem.MolFromSmarts(aldoxime_smarts)
    if pattern is None:
        return False, "Error in SMARTS pattern for aldoxime"
    
    # Search for the aldoxime substructure in the molecule.
    if mol.HasSubstructMatch(pattern):
        return True, "Contains an aldoxime group (RCH=NOH)"
    else:
        return False, "Aldoxime group not found in the molecule"

# Example usage:
# test_smiles = [
#     "[H]\\C(C(C)C)=N/O",  # (E)-2-methylpropanal oxime
#     "C([C@@H](/C(=N/O)/[H])C)C",  # (1E,2S)-2-methylbutanal oxime
#     "COc1cc(\\C=N/O)nc(c1)-c1ccccn1",  # (Z)-4-methoxy-2,2-bipyridine-6-carbaldehyde oxime
#     "COc1cc(\\C=N\\O)nc(-c2ccccn2)c1OC"  # caerulomycin C
# ]
#
# for smi in test_smiles:
#     result, reason = is_aldoxime(smi)
#     print(f"SMILES: {smi}\nResult: {result}\nReason: {reason}\n")