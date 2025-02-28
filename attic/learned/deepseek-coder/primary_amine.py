"""
Classifies: CHEBI:32877 primary amine
"""
"""
Classifies: CHEBI:32877 primary amine
"""
from rdkit import Chem

def is_primary_amine(smiles: str):
    """
    Determines if a molecule is a primary amine based on its SMILES string.
    A primary amine is a compound formally derived from ammonia by replacing one hydrogen atom by a hydrocarbyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a primary amine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for primary amine: -NH2 attached to any carbon
    primary_amine_pattern = Chem.MolFromSmarts("[NH2][CX4]")
    
    # Check if the molecule matches the primary amine pattern
    if mol.HasSubstructMatch(primary_amine_pattern):
        # Ensure that the nitrogen is not part of a secondary or tertiary amine
        for match in mol.GetSubstructMatches(primary_amine_pattern):
            nitrogen_atom = mol.GetAtomWithIdx(match[0])
            # Check that the nitrogen has exactly 2 hydrogens and 1 carbon neighbor
            if (nitrogen_atom.GetTotalNumHs() == 2 and 
                len([nbr for nbr in nitrogen_atom.GetNeighbors() if nbr.GetAtomicNum() == 6]) == 1):
                return True, "Contains a primary amine group (-NH2 attached to a carbon)"
        
    # Also check for aromatic amines (e.g., aniline derivatives)
    aromatic_amine_pattern = Chem.MolFromSmarts("[NH2][c]")
    if mol.HasSubstructMatch(aromatic_amine_pattern):
        for match in mol.GetSubstructMatches(aromatic_amine_pattern):
            nitrogen_atom = mol.GetAtomWithIdx(match[0])
            if nitrogen_atom.GetTotalNumHs() == 2:
                return True, "Contains a primary amine group (-NH2 attached to an aromatic carbon)"
    
    return False, "No primary amine group (-NH2 attached to a carbon) found"