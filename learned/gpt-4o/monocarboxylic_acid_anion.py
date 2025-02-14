"""
Classifies: CHEBI:35757 monocarboxylic acid anion
"""
from rdkit import Chem

def is_monocarboxylic_acid_anion(smiles: str):
    """
    Determines if a molecule is a monocarboxylic acid anion based on its SMILES string.
    A monocarboxylic acid anion has one carboxylate group, represented by the pattern C(=O)[O-],
    and lacks other charged functional groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monocarboxylic acid anion, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Identify carboxylate group pattern (C(=O)[O-])
    carboxylate_pattern = Chem.MolFromSmarts("C(=O)[O-]")
    carboxylate_matches = mol.GetSubstructMatches(carboxylate_pattern)
    
    # Check for exactly one carboxylate group
    if len(carboxylate_matches) != 1:
        if len(carboxylate_matches) == 0:
            return False, "No carboxylate group found"
        else:
            return False, f"Found {len(carboxylate_matches)} carboxylate groups, expected exactly 1"
    
    # Exclude molecules with additional anionic groups that may indicate more than monocarboxylate
    # This could be generalized as [O-] not part of carboxylate or other complex patterns
    non_carboxylate_anion_pattern = Chem.MolFromSmarts("[O-]~[!C]=[!O]")
    if mol.HasSubstructMatch(non_carboxylate_anion_pattern):
        return False, "Contains non-carboxylate anionic group, indicating it's not a simple monocarboxylic acid anion"

    return True, "Contains exactly one carboxylate group, indicating a monocarboxylic acid anion"