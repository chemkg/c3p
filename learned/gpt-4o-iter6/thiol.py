"""
Classifies: CHEBI:29256 thiol
"""
from rdkit import Chem

def is_thiol(smiles: str):
    """
    Determines if a molecule is a thiol based on its SMILES string.
    A thiol is defined as having a thiol group (-SH) attached to a carbon atom in aliphatic or aromatic moiety.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a thiol, False otherwise
        str: Reason for classification
    """
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the thiol group pattern (S-H)
    thiol_pattern = Chem.MolFromSmarts("[CX4,CX3,CX2,CX1]S")
    
    # Check for presence of thiol group
    if mol.HasSubstructMatch(thiol_pattern):
        return True, "Contains thiol group (-SH) attached to carbon"
    
    return False, "No thiol group attached to carbon found"