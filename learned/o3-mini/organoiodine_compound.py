"""
Classifies: CHEBI:37142 organoiodine compound
"""
"""
Classifies: Organoiodine compounds
Definition: An organoiodine compound is a compound containing at least one carbon-iodine bond.
"""

from rdkit import Chem

def is_organoiodine_compound(smiles: str):
    """
    Determines if a molecule is an organoiodine compound based on its SMILES string.
    An organoiodine compound has at least one carbon-iodine bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an organoiodine compound, False otherwise.
        str: Reason for classification.
    """
    # Parse the SMILES string into an RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern that matches a carbon-iodine bond
    # [#6] matches any carbon, and [I] matches iodine.
    pattern = Chem.MolFromSmarts("[#6]-[I]")
    
    # Check if the molecule contains at least one carbon-iodine bond
    if mol.HasSubstructMatch(pattern):
        return True, "Molecule contains at least one carbon-iodine bond"
    else:
        return False, "No carbon-iodine bond found in the molecule"