"""
Classifies: CHEBI:23003 carbamate ester
"""
"""
Classifies: carbamate ester
"""
from rdkit import Chem

def is_carbamate_ester(smiles: str):
    """
    Determines if a molecule is a carbamate ester based on its SMILES string.
    A carbamate ester has the general structure -[N]-C(=O)-O-, where nitrogen may have substituents.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a carbamate ester, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # SMARTS pattern for carbamate ester group
    carbamate_pattern = Chem.MolFromSmarts("[NX3;H0,H1,H2]-C(=O)-[OX2]")
    
    # Check for the carbamate pattern
    if mol.HasSubstructMatch(carbamate_pattern):
        return True, "Contains a carbamate ester group"
    else:
        return False, "Does not contain a carbamate ester group"