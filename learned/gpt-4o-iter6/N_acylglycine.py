"""
Classifies: CHEBI:16180 N-acylglycine
"""
from rdkit import Chem

def is_N_acylglycine(smiles: str):
    """
    Determines if a molecule is an N-acylglycine based on its SMILES string.
    An N-acylglycine is characterized by an acyl group attached to the nitrogen atom of glycine.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-acylglycine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for N-acylglycine pattern
    n_acylglycine_pattern = Chem.MolFromSmarts("[NX3][CX3](=[OX1])[CX3][OX2H]")
    if not mol.HasSubstructMatch(n_acylglycine_pattern):
        return False, "No N-acylglycine structure found"

    return True, "Molecule contains an N-acylglycine structure"