"""
Classifies: CHEBI:23437 cyanohydrin
"""
"""
Classifies: Cyanohydrin (alpha-hydroxynitrile)
Definition: An alpha-hydroxynitrile resulting from the formal addition of hydrogen cyanide to 
the C=O bond of an aldehyde or ketone. This program checks for a substructure where a tetrahedral 
carbon bears a free hydroxyl (–OH) group and is directly bonded (by a single bond) to a nitrile (C≡N) group.
"""

from rdkit import Chem

def is_cyanohydrin(smiles: str):
    """
    Determines if a molecule is a cyanohydrin based on its SMILES string.
    A cyanohydrin contains an alpha-hydroxynitrile group where one sp3 carbon is connected to a free -OH group 
    and to a nitrile group (-C≡N) through a single bond.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if the molecule is classified as a cyanohydrin, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a SMARTS pattern that requires:
    # 1. A tetrahedral carbon ([CX4]) with an attached free hydroxyl group ([OX2H]).
    # 2. This same carbon is single-bonded (-) to a nitrile carbon ([C]#[N]).
    # This pattern should not match cases where the oxygen is part of an ether or glycoside.
    cyanohydrin_pattern = Chem.MolFromSmarts("[CX4]([OX2H])-[C]#[N]")
    if cyanohydrin_pattern is None:
        return False, "SMARTS pattern for cyanohydrin could not be compiled"

    # Check if the molecule matches the cyanohydrin pattern
    if mol.HasSubstructMatch(cyanohydrin_pattern):
        return True, "Molecule contains a free alpha-hydroxynitrile (cyanohydrin) functional group"
    else:
        return False, "Alpha-hydroxynitrile (cyanohydrin) motif not detected in the molecule"