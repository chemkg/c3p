"""
Classifies: CHEBI:16366 anthocyanidin cation
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_anthocyanidin_cation(smiles: str):
    """
    Determines if a molecule is an anthocyanidin cation based on its SMILES string.
    An anthocyanidin cation is an aglycon of anthocyanin cation, characterized by a flavylium core.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an anthocyanidin cation, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the flavylium core SMARTS pattern, allowing for substitutions
    # The pattern is [O+]=C1c2ccccc2Cc3cc(O)cc(O)c13 to allow for substitutions,
    # This is a 2-phenylchromenylium ion with two OH substitutions.
    # This is a more flexible pattern to capture the core structure of anthocyanidins
    flavylium_core_pattern = Chem.MolFromSmarts("[O+]=C1c2ccccc2Cc3cc(O)cc(O)c13")


    # Check for the flavylium core
    if not mol.HasSubstructMatch(flavylium_core_pattern):
        return False, "Flavylium core not found"

    return True, "Contains the flavylium cation core (aglycone)"