"""
Classifies: CHEBI:27325 xanthophyll
"""
"""
Classifies: CHEBI:37208 xanthophyll

A subclass of carotenoids consisting of the oxygenated carotenes.
"""

from rdkit import Chem

def is_xanthophyll(smiles: str):
    """
    Determines if a molecule is a xanthophyll based on its SMILES string.
    A xanthophyll is an oxygenated carotenoid, containing at least one oxygen atom.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a xanthophyll, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for oxygen atoms
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if oxygen_count == 0:
        return False, "No oxygen atoms found, not a xanthophyll"
    
    # Check for conjugated carbon chain
    carbon_chain_pattern = Chem.MolFromSmarts("[cH2]=[c;!r]")
    carbon_chain_matches = mol.GetSubstructMatches(carbon_chain_pattern)
    if not carbon_chain_matches:
        return False, "No conjugated carbon chain found, not a carotenoid"
    
    return True, "Contains oxygenated conjugated carbon chain, characteristic of xanthophylls"