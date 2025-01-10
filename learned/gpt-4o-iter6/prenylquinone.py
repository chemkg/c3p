"""
Classifies: CHEBI:26255 prenylquinone
"""
from rdkit import Chem

def is_prenylquinone(smiles: str):
    """
    Determines if a molecule is a prenylquinone based on its SMILES string.
    A prenylquinone is a quinone substituted by a polyprenyl-derived side-chain.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a prenylquinone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define enhanced quinone patterns
    quinone_patterns = [ 
        Chem.MolFromSmarts("C1=CC(=O)C=C(C=O)C1"), # Benzoquinone
        Chem.MolFromSmarts("C1=CC2=C(C=C1)C(=O)C=CC2=O"), # Naphthoquinone
        Chem.MolFromSmarts("C1=CC(=O)C2=CC=CC(=O)C2=C1"), # Anthraquinone
        Chem.MolFromSmarts("O=C1C=CC(=O)C=C1") # Generic quinone
    ]
    
    # Check for presence of quinone core structure
    quinone_found = any(mol.HasSubstructMatch(pattern) for pattern in quinone_patterns)
    if not quinone_found:
        return False, "No quinone core structure found"

    # Define flexible prenyl chain pattern
    prenyl_chain_pattern = Chem.MolFromSmarts("C=C-C-C") # Basic isoprene unit
    # Consider branching and longer aliphatic prenyl chains
    extended_prenyl_pattern = Chem.MolFromSmarts("C(=C)CC")  # More flexible pattern for prenyl chains

    # Check for presence of prenyl or extended isoprene chains
    prenyl_matches = mol.GetSubstructMatches(prenyl_chain_pattern) + mol.GetSubstructMatches(extended_prenyl_pattern)
    if len(prenyl_matches) < 2:
        return False, "Insufficient prenyl-derived side chains found"
    
    return True, "Contains quinone structure with prenyl-derived side chain"