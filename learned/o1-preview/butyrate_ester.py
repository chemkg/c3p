"""
Classifies: CHEBI:50477 butyrate ester
"""
"""
Classifies: butyrate ester

Definition: Any carboxylic ester where the carboxylic acid component is butyric acid.
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_butyrate_ester(smiles: str):
    """
    Determines if a molecule is a butyrate ester based on its SMILES string.
    A butyrate ester is an ester where the carboxylic acid component is butyric acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a butyrate ester, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define butyrate ester pattern
    butyrate_ester_pattern = Chem.MolFromSmarts('[CH3][CH2][CH2]C(=O)O[*]')
    if mol.HasSubstructMatch(butyrate_ester_pattern):
        return True, "Contains butyrate ester group"

    return False, "Does not contain butyrate ester group"