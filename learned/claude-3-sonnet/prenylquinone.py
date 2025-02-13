"""
Classifies: CHEBI:26255 prenylquinone
"""
"""
Classifies: CHEBI:51714 prenylquinone
"""
from rdkit import Chem
from rdkit.Chem import rdqueries

def is_prenylquinone(smiles: str):
    """
    Determines if a molecule is a prenylquinone based on its SMILES string.
    A prenylquinone is a quinone substituted by a polyprenyl-derived side-chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a prenylquinone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for quinone substructure
    quinone_pattern = rdqueries.IsotonicQueriesFromSmarts("O=C1C(=O)C=CC=C1")
    if not mol.HasSubstructMatch(quinone_pattern):
        return False, "No quinone substructure found"

    # Check for prenyl side-chain(s)
    prenyl_pattern = rdqueries.IsotonicQueriesFromSmarts("[CH2]=[CH][CH2][CH2][CH2][CH3]")
    prenyl_matches = mol.GetSubstructMatches(prenyl_pattern)
    if not prenyl_matches:
        return False, "No prenyl side-chain found"

    # Check that prenyl chain(s) are connected to or near the quinone
    quinone_atoms = set(atom.GetIdx() for match in mol.GetSubstructMatches(quinone_pattern) for atom in match)
    prenyl_atoms = set(atom for match in prenyl_matches for atom in match)
    nearby_atoms = set()
    for prenyl_atom in prenyl_atoms:
        nearby_atoms.update(mol.GetAtomWithIdx(prenyl_atom).GetNeighbors())
    if not quinone_atoms.intersection(nearby_atoms):
        return False, "Prenyl side-chain(s) not connected or near the quinone"

    return True, "Contains a quinone substructure with one or more prenyl side-chains"