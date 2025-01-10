"""
Classifies: CHEBI:35746 fatty aldehyde
"""
"""
Classifies: CHEBI:63850 fatty aldehyde
"""
from rdkit import Chem

def is_fatty_aldehyde(smiles: str):
    """
    Determines if a molecule is a fatty aldehyde based on its SMILES string.
    A fatty aldehyde is defined as an aldehyde with a carbonyl group at one end
    of a carbon chain, typically derived from the reduction of a fatty acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty aldehyde, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for terminal aldehyde group pattern (C=O at terminal position)
    aldehyde_pattern = Chem.MolFromSmarts("[CX3H1](=O)[#6]")
    if not mol.HasSubstructMatch(aldehyde_pattern):
        return False, "No terminal aldehyde group found"

    # Check for carbon count, allowing for smaller aldehydes
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)

    # Allow smaller sizes, focusing on typical fatty chains
    if c_count < 4:
        return False, f"Too few carbons for a typical fatty aldehyde (found {c_count})"

    # Flexible pattern to count potential long chains, including unsaturation
    # Allow substructure like CCC=C, and permit branching
    flexible_chain_pattern = Chem.MolFromSmarts("[C&R1]~[C&R1]~[C&R1]~[C&R1]")
    if not mol.HasSubstructMatch(flexible_chain_pattern):
        return False, "Does not have a sufficient long carbon chain, considering unsaturations and branches"

    return True, "Contains terminal aldehyde group with a carbon chain suitable as fatty aldehyde"