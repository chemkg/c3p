"""
Classifies: CHEBI:33566 catechols
"""
from rdkit import Chem

def is_catechols(smiles: str):
    """
    Determines if a molecule is a catechol based on its SMILES string.
    A catechol is defined as a compound containing an o-diphenol (two hydroxyl groups on adjacent carbons) component.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a catechol, False otherwise
        str: Reason for classification
    """
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for o-diphenol: two hydroxyls on adjacent sp2 or sp3 carbons
    # This is the most critical improvement
    catechol_pattern = Chem.MolFromSmarts("[CX3,CX4]([OX2H,OX1-])[CX3,CX4]([OX2H,OX1-])")
    if catechol_pattern is None:
        return None, "Invalid SMARTS pattern"
    
    # Define SMARTS pattern for excluding tropolones
    tropolone_pattern = Chem.MolFromSmarts("O=C1C(O)C(O)=C[C,c][C,c]1")
    if tropolone_pattern is None:
        return None, "Invalid SMARTS pattern"

    # Check for substructure match for o-diphenol
    if mol.HasSubstructMatch(catechol_pattern):
        # check for the substructure match for tropolones
        if not mol.HasSubstructMatch(tropolone_pattern):
           return True, "Molecule contains an ortho-diphenol moiety"
        else:
            return False, "Molecule contains a tropolone structure, therefore not a catechol."    
    else:
        return False, "Molecule does not contain an ortho-diphenol moiety"