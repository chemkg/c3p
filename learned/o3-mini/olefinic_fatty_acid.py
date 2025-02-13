"""
Classifies: CHEBI:53339 olefinic fatty acid
"""
"""
Classifies: Olefinic Fatty Acid
Definition: Any fatty acid (i.e. molecule with a carboxylic acid group) that contains at least one C=C double bond.
"""

from rdkit import Chem

def is_olefinic_fatty_acid(smiles: str):
    """
    Determines if a molecule is an olefinic fatty acid based on its SMILES string.
    Olefinic fatty acids are defined as fatty acids (molecules containing a carboxylic acid group)
    that have at least one carbon–carbon double bond.
    
    Args:
        smiles (str): SMILES string of the molecule, e.g. "CCCC(=O)O" (fatty acid) with additional C=C bonds.
    
    Returns:
        bool: True if molecule is an olefinic fatty acid, False otherwise.
        str: Explanation of the classification decision.
    """
    # Parse the input SMILES string to a molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of a carboxylic acid group.
    # The SMARTS "C(=O)[O]" should capture carboxyl groups including those that may be represented in standard form.
    acid_pattern = Chem.MolFromSmarts("C(=O)[O]")
    if not mol.HasSubstructMatch(acid_pattern):
        return False, "No carboxylic acid group found; not a fatty acid"

    # Check for the presence of at least one carbon–carbon double bond.
    # The SMARTS "C=C" matches a carbon–carbon double bond (and avoids confusing it with a carbonyl C=O).
    olefin_pattern = Chem.MolFromSmarts("C=C")
    if not mol.HasSubstructMatch(olefin_pattern):
        return False, "No C=C double bonds found; not an olefinic fatty acid"

    return True, "Contains a carboxylic acid group and at least one C=C double bond, classifying it as an olefinic fatty acid"