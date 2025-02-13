"""
Classifies: CHEBI:50477 butyrate ester
"""
"""
Classifies: Butyrate ester (Any carboxylic ester where the carboxylic acid component is butyric acid)
Butyric acid is CH3CH2CH2COOH, so the corresponding ester moiety (after removal of OH) is CH3CH2CH2C(=O)O.
We use a SMARTS pattern that matches a linear chain: [CH3][CH2][CH2]C(=O)O.
"""

from rdkit import Chem

def is_butyrate_ester(smiles: str):
    """
    Determines if a molecule (given as a SMILES string) is a butyrate ester.
    A butyrate ester is defined as a carboxylic ester where the carboxylic acid component is butyric acid
    (i.e. the acyl group is CH3CH2CH2C(=O)-).

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule contains a butyrate ester fragment, False otherwise.
        str: Reason for classification.
    """
    
    # Parse the SMILES string into an RDKit molecule object.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS for an ester group in general.
    # But we only want the ester where the acyl (acid) component is n-butyric acid: CH3CH2CH2C(=O)O.
    # The SMARTS below matches a terminal methyl group attached to two CH2 groups then a carbonyl and an oxygen:
    butyrate_pattern = Chem.MolFromSmarts("[CH3][CH2][CH2]C(=O)O")
    if butyrate_pattern is None:
        return False, "Error in defining butyrate SMARTS pattern"
    
    # Check if the molecule has a substructure matching the butyrate ester fragment.
    if mol.HasSubstructMatch(butyrate_pattern):
        return True, "Molecule contains an ester group with the n-butyric acid moiety (CH3CH2CH2C(=O)O)."
    else:
        return False, "No butyrate ester substructure (CH3CH2CH2C(=O)O) found."
        
# Example usage (you can test each SMILES to see if it is classified as a butyrate ester)
if __name__ == "__main__":
    test_smiles = [
        "CCCC(=O)OCC",         # ethyl butyrate: should return True.
        "CCCC(=O)OC(C)CC",      # sec-butyl butyrate: should return True.
        "CC(=O)OC",            # methyl acetate (not butyrate): should return False.
    ]
    for sm in test_smiles:
        result, reason = is_butyrate_ester(sm)
        print(f"SMILES: {sm}")
        print(f"Result: {result}, Reason: {reason}\n")