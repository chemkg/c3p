"""
Classifies: CHEBI:17478 aldehyde
"""
"""
Classifies: Aldehyde (RC(=O)H)
Definition: A compound in which a carbonyl group (C=O) is bonded to one hydrogen and one R-group.
"""

from rdkit import Chem

def is_aldehyde(smiles: str):
    """
    Determines if a molecule is an aldehyde based on its SMILES string.
    An aldehyde, by definition, has the fragment RC(=O)H (a carbonyl group with exactly one hydrogen attached).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule contains an aldehyde functional group, False otherwise.
        str: Explanation for the classification.
    """
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # Define a SMARTS pattern for an aldehyde functional group.
    # The pattern "[CX3H1](=O)[!#1]" is interpreted as:
    #   [CX3H1] -> a trivalent (sp2 hybridized) carbon atom with exactly one hydrogen,
    #   (=O)    -> doubly bonded to an oxygen atom,
    #   [!#1]   -> and connected to an atom that is not hydrogen (the R group).
    aldehyde_smarts = "[CX3H1](=O)[!#1]"
    aldehyde_pattern = Chem.MolFromSmarts(aldehyde_smarts)
    
    if mol.HasSubstructMatch(aldehyde_pattern):
        return True, "Aldehyde group detected: contains a carbonyl (C=O) bonded to one hydrogen and one R-group."
    
    return False, "No aldehyde functional group (RC(=O)H) detected."

# Example of how one might test the function (uncomment to run):
# test_smiles = "OC=1C(=CC=CC1)/C=C/C=O"  # (E)-3-(2-Hydroxyphenyl)-2-propenal
# result, reason = is_aldehyde(test_smiles)
# print(result, reason)