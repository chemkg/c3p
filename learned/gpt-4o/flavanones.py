"""
Classifies: CHEBI:28863 flavanones
"""
from rdkit import Chem

def is_flavanones(smiles: str):
    """
    Determines if a molecule is a flavanone based on its SMILES string.
    A flavanone is defined by a 3,4-dihydro-2-aryl-2H-1-benzopyran-4-one skeleton,
    and its substituted derivatives.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a flavanone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for the flavanone core
    # This pattern is for a 3,4-dihydro-2-aryl-2H-1-benzopyran-4-one
    flavanone_pattern = Chem.MolFromSmarts("O=C1CC2=C(O1)C=CC=C2c3ccccc3")

    if not mol.HasSubstructMatch(flaНавanone_pattern):
        return False, "No 3,4-dihydro-2-aryl-2H-1-benzopyran-4-one skeleton found"

    return True, "Contains the 3,4-dihydro-2-aryl-2H-1-benzopyran-4-one skeleton"