"""
Classifies: CHEBI:37142 organoiodine compound
"""
from rdkit import Chem

def is_organoiodine_compound(smiles: str):
    """
    Determines if a molecule is an organoiodine compound based on its SMILES string.
    An organoiodine compound should contain at least one carbon-iodine bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an organoiodine compound, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carbon-iodine bond pattern using SMARTS
    c_i_pattern = Chem.MolFromSmarts("[#6]-[#53]")

    # Split the SMILES string for potential disconnected components and check each
    fragments = Chem.GetMolFrags(mol, asMols=True)
    for fragment in fragments:
        if fragment.HasSubstructMatch(c_i_pattern):
            return True, "Contains carbon-iodine bond"

    return False, "No carbon-iodine bond found"