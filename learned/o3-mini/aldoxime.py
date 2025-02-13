"""
Classifies: CHEBI:22307 aldoxime
"""
"""
Classifies: Aldoxime (Oximes of aldehydes RCH=NOH)

An aldoxime is an oxime derivative of an aldehyde, meaning that the oxime group is attached 
to a carbon that originally was part of an aldehyde (RCH=O) and should still display a single hydrogen on it.
This script checks for a substructure match of a CH group double-bonded to a nitrogen that is bound to an OH.
"""

from rdkit import Chem

def is_aldoxime(smiles: str):
    """
    Determines if a molecule is an aldoxime based on its SMILES string.
    An aldoxime is an oxime derived from an aldehyde and has the functional group RCH=NOH.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is classified as an aldoxime, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string into an RDKit molecule object.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Add explicit hydrogens so we can check for the presence of exactly one hydrogen on the carbon.
    mol = Chem.AddHs(mol)

    # Define a SMARTS pattern for an aldoxime.
    # Pattern explanation:
    # [#6H1] - a carbon atom that has exactly one hydrogen (typical of an aldehyde carbon).
    # =      - double bond.
    # [NX2]  - a nitrogen atom with valence 2 (in the C=N bond).
    # [OX1H] - an oxygen atom with one hydrogen (the hydroxyl group).
    pattern = Chem.MolFromSmarts("[#6H1]=[NX2][OX1H]")
    if pattern is None:
        return False, "Error in SMARTS pattern"

    # Check if the molecule contains the aldoxime substructure.
    if mol.HasSubstructMatch(pattern):
        return True, "Contains an aldoxime group (RCH=NOH)"
    else:
        return False, "Aldoxime group not found in the molecule"