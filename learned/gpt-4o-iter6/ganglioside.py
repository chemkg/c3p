"""
Classifies: CHEBI:28892 ganglioside
"""
from rdkit import Chem

def is_ganglioside(smiles: str):
    """
    Determines if a molecule is a ganglioside based on its SMILES string.
    A ganglioside is composed of a glycosphingolipid with one or more sialic acids.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a ganglioside, False otherwise
        str: Reason for classification
    """
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Generalized ceramide backbone pattern 
    ceramide_pattern = Chem.MolFromSmarts("O[C@@H](CO)C[NH]C(=O)C(CCCCC)O[C@H](CCCC)")

    # Generalized sialic acid pattern allowing for N-acetyl groups and other variations
    sialic_acid_pattern = Chem.MolFromSmarts("C[C@@H](O)[C@H](C=O)OC(CO)COC(=O)")

    # Check for generalized ceramide backbone presence
    if not mol.HasSubstructMatch(ceramide_pattern):
        return False, "No generalized ceramide backbone found"

    # Check for at least one sialic acid residue
    sialic_acid_residues = mol.GetSubstructMatches(sialic_acid_pattern)
    if len(sialic_acid_residues) < 1:
        return False, "No sialic acid residues found"

    return True, "Contains generalized ceramide backbone and one or more sialic acids"