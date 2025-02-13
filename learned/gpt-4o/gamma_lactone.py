"""
Classifies: CHEBI:37581 gamma-lactone
"""
from rdkit import Chem

def is_gamma_lactone(smiles: str):
    """
    Determines if a molecule is a gamma-lactone based on its SMILES string.
    A gamma-lactone is characterized by the presence of a five-membered lactone ring.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is a gamma-lactone, False otherwise.
        str: Explanation for the classification decision.
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a SMARTS pattern for a five-membered gamma-lactone ring.
    # Note: This pattern includes the possibility of substituents, and allows for stereochemistry.
    gamma_lactone_pattern = Chem.MolFromSmarts("O=C1OC[CH2][CH2]1")

    if mol.HasSubstructMatch(gamma_lactone_pattern):
        return True, "Contains a five-membered lactone ring characteristic of gamma-lactones"
    else:
        return False, "No five-membered lactone ring found"