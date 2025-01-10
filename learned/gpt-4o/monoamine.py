"""
Classifies: CHEBI:63534 monoamine
"""
from rdkit import Chem

def is_monoamine(smiles: str):
    """
    Determines if a molecule is classified as a monoamine based on its SMILES string.
    A monoamine contains an amino group connected to an aromatic ring by a two-carbon chain.

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

    # Updated aromatic ring pattern to include heteroaromatics
    aromatic_ring_pattern = Chem.MolFromSmarts("a1aaaaa1")  # General pattern for any aromatic ring
    if not mol.HasSubstructMatch(aromatic_ring_pattern):
        return False, "No appropriate aromatic ring found"

    # Update amine pattern to capture all forms specifically NX3 ammonia-type
    amine_pattern = Chem.MolFromSmarts("[NX3,NX4+][C,c]")  # Ensure amine is connected to a carbon or aromatic carbon
    if not mol.HasSubstructMatch(amine_pattern):
        return False, "No amino group found"

    # Refined two-carbon chain linking amine to aromatic ring
    chain_pattern = Chem.MolFromSmarts("[a][CX4][CX4][NX3,NX4+]")  # Aromatic carbon -> aliphatic C -> C -> amine
    if not mol.HasSubstructMatch(chain_pattern):
        return False, "No two-carbon chain connecting amino to aromatic ring found"

    return True, "Contains amino group connected to aromatic ring by a two-carbon chain"

# Example debugging
# print(is_monoamine("CNC[C@H](O)c1ccc(O)c(O)c1"))  # Expected: True
# print(is_monoamine("NCCc1ccc(O)c(O)c1"))  # Expected: True