"""
Classifies: CHEBI:38958 indole alkaloid
"""
"""
Classifies: CHEBI:33831 indole alkaloid
"An alkaloid containing an indole skeleton."
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_indole_alkaloid(smiles: str):
    """
    Determines if a molecule is an indole alkaloid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an indole alkaloid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for indole substructure that may be part of a larger ring system
    indole_alkaloid_pattern = Chem.MolFromSmarts("c1ccc2c(c1)cnc2")
    has_indole_alkaloid = mol.HasSubstructMatch(indole_alkaloid_pattern)

    # Look for basic nitrogen
    basic_nitrogen_pattern = Chem.MolFromSmarts("[N;H0;+]")
    has_basic_nitrogen = mol.HasSubstructMatch(basic_nitrogen_pattern)

    # Classify as indole alkaloid if both conditions are met
    if has_indole_alkaloid and has_basic_nitrogen:
        return True, "Contains an indole skeleton fused to a larger ring system and a basic nitrogen atom"
    elif not has_indole_alkaloid:
        return False, "Does not contain an indole skeleton fused to a larger ring system"
    else:
        return False, "Does not contain a basic nitrogen atom"