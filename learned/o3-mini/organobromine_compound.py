"""
Classifies: CHEBI:37141 organobromine compound
"""
"""
Classifies: Organobromine compound
Definition: A compound containing at least one carbon-bromine bond.
"""

from rdkit import Chem

def is_organobromine_compound(smiles: str):
    """
    Determines if a molecule is an organobromine compound based on its SMILES string.
    An organobromine compound contains at least one carbon-bromine bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an organobromine compound, False otherwise.
        str: Reason for classification.
    """
    # Parse the SMILES string to a molecule object.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a SMARTS pattern for a carbon-bromine bond.
    # [#6] matches any carbon, and Br matches a bromine atom.
    c_br_pattern = Chem.MolFromSmarts("[#6]-Br")
    
    # Check for at least one carbon-bromine bond in the molecule.
    if mol.HasSubstructMatch(c_br_pattern):
        return True, "Contains at least one carbon-bromine bond."
    else:
        return False, "No carbon-bromine bond found."