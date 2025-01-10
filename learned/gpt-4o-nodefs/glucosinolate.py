"""
Classifies: CHEBI:24279 glucosinolate
"""
from rdkit import Chem

def is_glucosinolate(smiles: str):
    """
    Determines if a molecule is a glucosinolate based on its SMILES string.
    Glucosinolates contain a part of a glucose-derived moiety, a sulfur linkage, 
    and an N-sulfooxyimino group.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glucosinolate, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns
    # Relaxed partial glucose pattern to capture more variation:
    glucose_pattern = Chem.MolFromSmarts("[C@H]1(O[C@@H]([C@H](O)[C@@H](O)[C@H](O)C1)CO)")

    # General sulfur linkage pattern, emphasizing flexible placement:
    sulfur_linkage_pattern = Chem.MolFromSmarts("[SX2]-*")

    # N-sulfooxyimino group pattern:
    sulfooxyimino_pattern = Chem.MolFromSmarts("N=OS(=O)(=O)[O-]")

    # Check for essential glucosinolate components
    if not mol.HasSubstructMatch(glucose_pattern):
        return False, "No recognizable part-glucose scaffold found"
    
    if not mol.HasSubstructMatch(sulfur_linkage_pattern):
        return False, "No identifiable sulfur linkage found"

    if not mol.HasSubstructMatch(sulfooxyimino_pattern):
        return False, "Missing N-sulfooxyimino functional group"

    # If all key components are present, classify as glucosinolate
    return True, "Contains part-glucose scaffold, sulfur linkage, and N-sulfooxyimino group consistent with glucosinolate"