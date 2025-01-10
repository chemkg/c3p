"""
Classifies: CHEBI:73155 trienoic fatty acid
"""
from rdkit import Chem

def is_trienoic_fatty_acid(smiles: str):
    """
    Determines if a molecule is a trienoic fatty acid based on its SMILES string.
    A trienoic fatty acid is defined as a polyunsaturated fatty acid with exactly three double bonds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a trienoic fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group: It's a typical feature of fatty acids
    # but not the only form (e.g., esterified forms exist)
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[O]")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid or ester group found; not a simple fatty acid form"

    # Count all C=C double bonds generally, no need to be consecutive
    double_bond_pattern = Chem.MolFromSmarts("C=C") 
    double_bond_matches = len(mol.GetSubstructMatches(double_bond_pattern))

    if double_bond_matches == 3:
        return True, "Contains exactly three double bonds; classified as trienoic fatty acid"
    else:
        return False, f"Has {double_bond_matches} double bonds; not a trienoic fatty acid"