"""
Classifies: CHEBI:47916 flavonoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_flavonoid(smiles: str):
    """
    Determines if a molecule is a flavonoid based on its SMILES string.
    A flavonoid has a 1-benzopyran core with an aryl substituent at position 2.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a flavonoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS for the 1-benzopyran core with aryl substituent at position 2
    # This pattern captures a 6-membered ring fused to a 6-membered ring via an oxygen atom with a substituent at the position next to the O
    flavonoid_smarts = '[c]1[c](~[c]2[c]([c][c][c][c]2)[o][c]3[c][c][c][c][c]3)[c][c][c][c]1'


    # Convert SMARTS to Mol object
    flavonoid_mol = Chem.MolFromSmarts(flavonoid_smarts)

    # Check if the molecule matches the pattern
    if mol.HasSubstructMatch(flavonoid_mol):
        return True, "Has 1-benzopyran core with aryl substituent at position 2"
    else:
        return False, "No 1-benzopyran core with aryl substituent at position 2 found"