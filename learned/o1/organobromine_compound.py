"""
Classifies: CHEBI:37141 organobromine compound
"""
from rdkit import Chem

def is_organobromine_compound(smiles: str):
    """
    Determines if a molecule is an organobromine compound based on its SMILES string.
    An organobromine compound is defined as a compound containing at least one carbon-bromine bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an organobromine compound, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string into a molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a SMARTS pattern for a carbon-bromine bond
    carbon_bromine_pattern = Chem.MolFromSmarts("[#6]-Br")  # [#6] is carbon, Br is bromine

    # Check if the molecule contains a carbon-bromine bond
    if mol.HasSubstructMatch(carbon_bromine_pattern):
        return True, "Contains at least one carbon-bromine bond"
    else:
        return False, "Does not contain any carbon-bromine bonds"