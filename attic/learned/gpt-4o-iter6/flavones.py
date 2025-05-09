"""
Classifies: CHEBI:24043 flavones
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_flavones(smiles: str) -> (bool, str):
    """
    Determines if a molecule is a flavone based on its SMILES string.
    A flavone has a 2-aryl-1-benzopyran-4-one skeleton and its substituted derivatives.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a flavone, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a more flexible core structure for flavones: 2-aryl-1-benzopyran-4-one
    flavone_core_pattern = Chem.MolFromSmarts("c1cc(ccc1)-c2cc(=O)c3ccccc3o2")
    
    # Ensure pattern was defined correctly
    if flavone_core_pattern is None:
        return None, "Failed to construct flavone core pattern"

    # Ensure molecule has flavone core structure
    if not mol.HasSubstructMatch(flavone_core_pattern):
        return False, "No flavone core structure found"

    return True, "Contains flavone core structure"