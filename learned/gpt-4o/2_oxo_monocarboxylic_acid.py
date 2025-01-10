"""
Classifies: CHEBI:35910 2-oxo monocarboxylic acid
"""
from rdkit import Chem

def is_2_oxo_monocarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a 2-oxo monocarboxylic acid based on its SMILES string.
    A 2-oxo monocarboxylic acid has a 2-oxo substituent and a monocarboxylic acid group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2-oxo monocarboxylic acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS pattern for 2-oxo group: [CX3](=O)[#6] where #6 is a carbon
    oxo_pattern = Chem.MolFromSmarts("[#6][CX3](=O)[#6]")
    
    # SMARTS pattern for carboxylic acid: C(=O)O
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")

    # Check for 2-oxo group
    if not mol.HasSubstructMatch(oxo_pattern):
        return False, "No 2-oxo substituent found"

    # Check for carboxylic acid group
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # If both patterns are present, it's a 2-oxo monocarboxylic acid
    return True, "Contains both 2-oxo group and carboxylic acid group"