"""
Classifies: CHEBI:60834 beta-carbolines
"""
from rdkit import Chem

def is_beta_carbolines(smiles: str):
    """
    Determines if a molecule is a beta-carboline based on its SMILES string.
    Beta-carbolines typically contain a pyrido[3,4-b]indole or related structure.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a beta-carboline, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES to rdkit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a wider beta-carboline-like pattern potentially capturing variations
    beta_carboline_pattern = Chem.MolFromSmarts("c1ccc2c(c1)[nH]c3cccnc23")

    # Check if the molecule contains the beta-carboline core structure
    if mol.HasSubstructMatch(beta_carboline_pattern):
        return True, "Contains beta-carboline (9H-pyrido[3,4-b]indole or similar) core structure"
    else:
        return False, "Missing beta-carboline core structure"

# Example usage: 
# result, reason = is_beta_carbolines("CN1C2=C(C=CC(=C2)OC)C3=C1[C@H](N(CC34CCN(CC4)C(=O)C5=CC=CC=N5)C(=O)COC)CO")
# print(result, reason)