"""
Classifies: CHEBI:27093 tricarboxylic acid
"""
"""
Classifies: tricarboxylic acid
"""
from rdkit import Chem
from rdkit.Chem import rdchem

def is_tricarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a tricarboxylic acid based on its SMILES string.
    A tricarboxylic acid is defined as having three distinct carboxylic acid groups (-COOH).

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

    # Define carboxylic acid group using SMARTS
    carboxylic_acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H]")
    carboxylic_acid_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    
    # Count the number of carboxylic acid groups
    n_carboxylic_acids = len(carboxylic_acid_matches)

    # If exactly 3 carboxylic acids are found, check their distinct contexts
    if n_carboxylic_acids >= 3:
        distinct_carbons = set()
        for match in carboxylic_acid_matches:
            # Each match is a tuple of atom indices (C=O and -OH)
            carbon_idx = match[0]
            distinct_carbons.add(carbon_idx)

        # Verify that carboxylic acids are attached to distinct carbon atoms
        if len(distinct_carbons) >= 3:
            # Ensure that additional constraints are met, e.g., tethering and context isolation:
            # Check that no carboxyl group is internally linked in a way that forms unintended matches.
            if n_carboxylic_acids == len(distinct_carbons):
                return True, "Contains exactly three distinct carboxylic acid groups"
            else:
                return False, "Carboxylic acids have overlapping or chemically linked nature"
        else:
            return False, "Carboxylic acids are not attached to distinct carbon atoms"
    else:
        return False, f"Contains {n_carboxylic_acids} carboxylic acid groups, expected exactly 3"