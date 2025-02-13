"""
Classifies: CHEBI:87658 decanoate ester
"""
"""
Classifies: Decanoate ester
A decanoate ester is defined as a fatty acid ester resulting from the formal condensation 
of the carboxy group of decanoic (capric) acid with the hydroxy group of an alcohol or phenol.
Decanoic acid is CH3-(CH2)8-C(=O)OH and thus the ester moiety will contain the acyl fragment:
CH3-(CH2)8-C(=O)O (or its deprotonated form).
"""
from rdkit import Chem

def is_decanoate_ester(smiles: str):
    """
    Determines if a molecule is a decanoate ester based on its SMILES string.
    The method searches for a decanoate acyl fragment:
      CH3-(CH2)8-C(=O)O  (or with a negatively charged oxygen instead of O)
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the decanoate ester substructure is found, False otherwise.
        str: Explanation or reason for the classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns for decanoate ester acyl fragment.
    # The first pattern matches CH3-CH2-CH2-CH2-CH2-CH2-CH2-CH2-CH2-C(=O)O
    # (note: 1 CH3 followed by eight CH2 groups gives 9 carbons; the carbonyl carbon is the 10th).
    decanoate_pattern1 = Chem.MolFromSmarts("[CH3][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2]C(=O)O")
    # The second pattern allows for the ester oxygen being deprotonated: [O-]
    decanoate_pattern2 = Chem.MolFromSmarts("[CH3][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2]C(=O)[O-]")

    # Check if either pattern is present in the molecule
    if mol.HasSubstructMatch(decanoate_pattern1) or mol.HasSubstructMatch(decanoate_pattern2):
        return True, "Contains decanoate ester moiety (decanoic acyl group detected)"
    else:
        return False, "Decanoate ester fragment not found in the molecule"