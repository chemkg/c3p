"""
Classifies: CHEBI:134251 guaiacols
"""
"""
Classifies: Guaiacols
Defined as: Any phenol carrying an additional methoxy substituent at the ortho-position.
"""
from rdkit import Chem

def is_guaiacols(smiles: str):
    """
    Determines if a molecule is a guaiacol based on its SMILES string.
    A guaiacol is defined as any phenol molecule carrying an additional methoxy group (–OCH3) on an ortho-position relative to the –OH.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule contains a guaiacol moiety, False otherwise.
        str: Reason for classification.
    """
    # Parse the SMILES string into a molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define two SMARTS patterns to capture a guaiacol moiety:
    # Pattern 1: aromatic ring with -OH in one position and -OCH3 in an adjacent position.
    pattern1 = Chem.MolFromSmarts("c1c(O)ccc(OC)c1")
    # Pattern 2: the substituents reversed in order.
    pattern2 = Chem.MolFromSmarts("c1c(OC)ccc(O)c1")

    # Check if either substructure pattern is found in the molecule.
    if mol.HasSubstructMatch(pattern1) or mol.HasSubstructMatch(pattern2):
        return True, "Contains a guaiacol moiety (phenol with an ortho-methoxy substituent)"
    else:
        return False, "Does not contain a guaiacol moiety: missing phenol with ortho-methoxy substituent"