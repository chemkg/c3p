"""
Classifies: CHEBI:139588 alpha-hydroxy ketone
"""
"""
Classifies: Alpha-hydroxy ketone
Definition: A ketone containing a hydroxy group on the alpha-carbon relative to the C=O group.
"""

from rdkit import Chem

def is_alpha_hydroxy_ketone(smiles: str):
    """
    Determines if a molecule is an alpha-hydroxy ketone based on its SMILES string.
    An alpha-hydroxy ketone is a ketone having an -OH group on one of the carbons directly adjacent
    to the carbonyl (C=O) group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an alpha-hydroxy ketone, False otherwise.
        str: Explanation of the classification.
    """
    # Parse SMILES into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns to capture the alpha-hydroxy ketone motif.
    # Pattern 1: A carbon bearing an -OH group adjacent to a carbonyl carbon
    # followed by any carbon. This represents: [C(–OH)]–C(=O)–C
    pattern1 = Chem.MolFromSmarts("[#6]([OX2H])-[CX3](=O)-[#6]")
    # Pattern 2: A carbon on the left, connected to a carbonyl, whose right neighbor bears an -OH.
    # This represents: C–C(=O)–[C(–OH)]
    pattern2 = Chem.MolFromSmarts("[#6]-[CX3](=O)-[#6]([OX2H])")

    # Check if either pattern is present in the molecule.
    if mol.HasSubstructMatch(pattern1) or mol.HasSubstructMatch(pattern2):
        return True, "Molecule contains a ketone with a hydroxy group on the alpha-carbon (alpha-hydroxy ketone detected)."
    else:
        return False, "No alpha-hydroxy ketone functional group found (ketone lacking an adjacent -OH group)."