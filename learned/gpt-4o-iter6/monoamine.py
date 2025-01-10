"""
Classifies: CHEBI:63534 monoamine
"""
from rdkit import Chem

def is_monoamine(smiles: str):
    """
    Determines if a molecule is a monoamine based on its SMILES string.
    A monoamine is defined as an arylamino compound which contains one amino group connected
    to an aromatic ring by a two-carbon chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monoamine, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Updated SMARTS pattern for an amino group connected to an aromatic system via a two-carbon chain
    # The pattern: [NX3;H2] for primary amine, [CH2] for methylene, and [c] for aromatic carbon.
    two_carbon_amino_pattern = Chem.MolFromSmarts("[NX3;H2]-[CH2]-[CH2]-[c]")

    # Check for two-carbon chain with amino group connected to aromatic system
    if mol.HasSubstructMatch(two_carbon_amino_pattern):
        return True, "Contains amino group connected to aromatic ring by a two-carbon chain"

    return False, "Missing two-carbon chain with amino group connected to an aromatic ring"