"""
Classifies: CHEBI:17792 organohalogen compound
"""
from rdkit import Chem

def is_organohalogen_compound(smiles: str):
    """
    Determines if a molecule is an organohalogen compound based on its SMILES string.
    An organohalogen compound contains at least one carbon-halogen bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an organohalogen compound, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define halogen atomic numbers
    halogens = {9, 17, 35, 53}  # F, Cl, Br, I

    # Traverse atoms to find carbon-halogen bonds
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() in halogens:
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 6:  # Carbon
                    return True, "Contains carbon-halogen bond"

    return False, "No carbon-halogen bond found"