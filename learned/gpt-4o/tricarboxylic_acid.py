"""
Classifies: CHEBI:27093 tricarboxylic acid
"""
from rdkit import Chem

def is_tricarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a tricarboxylic acid based on its SMILES string.
    A tricarboxylic acid contains three carboxy groups (-C(=O)OH).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a tricarboxylic acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for carboxylic acid pattern (C(=O)O)
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    carboxy_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)

    # Check if there are exactly three carboxyl groups
    if len(carboxy_matches) == 3:
        return True, "Contains three carboxylic acid groups"
    else:
        return False, f"Found {len(carboxy_matches)} carboxylic acid groups, need exactly 3"

# Example usage
smiles = "N[C@@H](CCC[C@H](NC(=O)CCC(O)=O)C(O)=O)C(O)=O"
print(is_tricarboxylic_acid(smiles))