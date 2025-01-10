"""
Classifies: CHEBI:31488 N-acylsphinganine
"""
from rdkit import Chem

def is_N_acylsphinganine(smiles: str):
    """
    Determines if a molecule is an N-acylsphinganine based on its SMILES string,
    which is a ceramide consisting of sphinganine in which one of the amino hydrogens
    is substituted by a fatty acyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-acylsphinganine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return False, "Invalid SMILES string"

    # Look for the sphinganine backbone pattern (including specific stereochemistry)
    sphinganine_pattern = Chem.MolFromSmarts("[C@H](O)[C@@H](CO)N")
    if sphinganine_pattern is None:
        return False, "Failed to build sphinganine pattern"
        
    if not mol.HasSubstructMatch(sphinganine_pattern):
        return False, "No sphinganine backbone found"

    # Look for N-acyl amide bond pattern (-NC(=O)-)
    amide_pattern = Chem.MolFromSmarts("NC(=O)")
    if amide_pattern is None:
        return False, "Failed to build amide pattern"
        
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "No amide bond (N-acyl group) found"

    # Check for a long aliphatic carbon chain (8 or more carbons) attached to the amide
    fatty_acid_chain_pattern = Chem.MolFromSmarts("C{8,}")
    if fatty_acid_chain_pattern is None:
        return False, "Failed to build fatty acid chain pattern"
        
    if not mol.HasSubstructMatch(fatty_acid_chain_pattern):
        return False, "No long aliphatic carbon chain found"

    return True, "Molecule is an N-acylsphinganine with sphinganine backbone and N-acyl group"