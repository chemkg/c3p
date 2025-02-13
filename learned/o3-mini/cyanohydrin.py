"""
Classifies: CHEBI:23437 cyanohydrin
"""
"""
Classifies: Cyanohydrin (alpha-hydroxynitrile)
Definition: An alpha-hydroxynitrile resulting from the formal addition of hydrogen cyanide to the C=O bond of an aldehyde or ketone.
This program checks for a substructure where a tetrahedral carbon bears a hydroxyl (–OH) group and is directly bonded to a nitrile (–C≡N) group.
"""

from rdkit import Chem

def is_cyanohydrin(smiles: str):
    """
    Determines if a molecule is a cyanohydrin based on its SMILES string.
    A cyanohydrin contains an alpha-hydroxynitrile group where one carbon has both a hydroxyl (OH) and a nitrile (C#N) substituent.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a cyanohydrin, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a SMARTS pattern that represents the alpha-hydroxynitrile substructure.
    # The pattern looks for a tetrahedral carbon ([CX4]) attached to an -OH group (O)
    # and bonded to a carbon that is part of a nitrile (C#N).
    cyanohydrin_pattern = Chem.MolFromSmarts("[CX4](O)[C]#[N]")
    if cyanohydrin_pattern is None:
        return False, "SMARTS pattern for cyanohydrin could not be compiled"

    # Check if the molecule has at least one match for the cyanohydrin pattern.
    if mol.HasSubstructMatch(cyanohydrin_pattern):
        return True, "Molecule contains an alpha-hydroxynitrile (cyanohydrin) functional group"
    else:
        return False, "Alpha-hydroxynitrile (cyanohydrin) motif not detected in the molecule"