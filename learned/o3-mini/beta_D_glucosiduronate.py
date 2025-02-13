"""
Classifies: CHEBI:83411 beta-D-glucosiduronate
"""
"""
Classifies: beta-D-glucosiduronate
A carbohydrate acid derivative anion obtained by deprotonation of the carboxy group of any beta-D-glucuronic acid.
Major species at pH 7.3.
"""

from rdkit import Chem

def is_beta_D_glucosiduronate(smiles: str):
    """
    Determines if a molecule is a beta-D-glucosiduronate based on its SMILES string.
    Beta-D-glucosiduronate is defined as the anion derived from beta-D-glucuronic acid.
    This function looks for a pyranose ring with defined beta stereochemistry and a carboxylate group,
    i.e. a C(=O)[O-] substituent at the appropriate position.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule contains a beta-D-glucuronide moiety, False otherwise.
        str: A reason for the classification.
    """
    # Parse the input SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a SMARTS pattern for the beta-D-glucuronide moiety.
    # This SMARTS looks for a pyranose ring with a carboxylate group.
    # Note: The stereochemistry is specified in the SMARTS to help distinguish the beta configuration.
    glucuronide_smarts = "O[C@@H]1[C@H](O)[C@@H](O)[C@@H](O)[C@H](O1)C(=O)[O-]"
    glucuronide_pattern = Chem.MolFromSmarts(glucuronide_smarts)
    if glucuronide_pattern is None:
        return False, "Error in SMARTS pattern definition"

    # Check if the molecule has the beta-D-glucuronide substructure.
    if not mol.HasSubstructMatch(glucuronide_pattern):
        return False, "Beta-D-glucuronide moiety not found in the molecule"

    # If desired, additional checks like overall charge or other criteria might be applied.
    return True, "Contains beta-D-glucuronide moiety (anionic form of beta-D-glucuronic acid)"

# For testing purposes (uncomment the lines below to test locally):
# test_smiles = "C1[C@@]2([C@]3(CC[C@]4([C@]([C@@]3([C@@H](C[C@@]2(C[C@@H](C1)O)[H])O)[H])(CC[C@@]4([C@@H](CCC(O[C@@H]5O[C@@H]([C@H]([C@@H]([C@H]5O)O)O)C([O-])=O)=O)C)[H])[H])C)[H])C"
# result, reason = is_beta_D_glucosiduronate(test_smiles)
# print(result, reason)