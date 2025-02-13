"""
Classifies: CHEBI:87691 tetradecanoate ester
"""
"""
Classifies: tetradecanoate ester
A tetradecanoate ester is defined as a fatty acid ester obtained by condensation
of the carboxy group of tetradecanoic acid (myristic acid) with the hydroxy group
of an alcohol or phenol.
This program checks if the given SMILES contains an ester substructure featuring the
tetradecanoate (myristate) moiety, i.e. a straight chain acyl group with 14 carbons (including
the carbonyl carbon) attached by an oxygen as part of an ester function.
"""

from rdkit import Chem

def is_tetradecanoate_ester(smiles: str):
    """
    Determines if a molecule is a tetradecanoate ester based on its SMILES string.
    
    The detection is based on the presence of an ester fragment derived from tetradecanoic acid,
    i.e. a group matching "CCCCCCCCCCCCCC(=O)O". This pattern corresponds to CH3-(CH2)12-COO-.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule contains a tetradecanoate ester moiety, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern for a tetradecanoate ester group.
    # "CCCCCCCCCCCCCC(=O)O" represents a 14-carbon chain (CH3-(CH2)12-CO) linked via an oxygen.
    myristate_pattern = Chem.MolFromSmarts("CCCCCCCCCCCCCC(=O)O")
    if myristate_pattern is None:
        return False, "Error in defining the tetradecanoate ester pattern"

    # Check if the molecule has a matching substructure.
    if mol.HasSubstructMatch(myristate_pattern):
        return True, "Contains a tetradecanoate ester moiety (myristate ester group detected)"
    else:
        return False, "No tetradecanoate ester moiety detected"
        
# Example usage (optional):
if __name__ == "__main__":
    # Test with one of the provided examples, e.g. "tetradecanoyl tetradecanoate"
    test_smiles = "C(CCCCCCCC)CCCCC(OCCCCCCCCCCCCCC)=O"
    result, reason = is_tetradecanoate_ester(test_smiles)
    print(result, reason)