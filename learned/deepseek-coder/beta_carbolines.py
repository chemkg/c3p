"""
Classifies: CHEBI:60834 beta-carbolines
"""
"""
Classifies: CHEBI:38830 beta-carboline
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_beta_carbolines(smiles: str):
    """
    Determines if a molecule is a beta-carboline based on its SMILES string.
    A beta-carboline is a pyrido[3,4-b]indole skeleton, which can be hydrogenated or substituted.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a beta-carboline, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the beta-carboline core structure (pyrido[3,4-b]indole)
    # This pattern is more flexible and accounts for substitutions and hydrogenations
    beta_carboline_pattern = Chem.MolFromSmarts("[nX2]1c2ccccc2c3c1cc[nX2,cX3]3")
    
    # Check if the molecule contains the beta-carboline core
    if mol.HasSubstructMatch(beta_carboline_pattern):
        return True, "Contains the beta-carboline core structure (pyrido[3,4-b]indole)"
    
    # If the core is not found, check for hydrogenated derivatives
    # This pattern accounts for fully or partially hydrogenated rings
    hydrogenated_pattern = Chem.MolFromSmarts("[nX2]1c2cccc[cX3,cX2]2c3c1cc[nX2,cX3]3")
    if mol.HasSubstructMatch(hydrogenated_pattern):
        return True, "Contains a hydrogenated beta-carboline core structure"
    
    # If no match is found, return False
    return False, "No beta-carboline core structure found"