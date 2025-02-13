"""
Classifies: CHEBI:47622 acetate ester
"""
"""
Classifies: Any carboxylic ester where the carboxylic acid component is acetic acid.
The function is_acetate_ester determines if a molecule (given by a SMILES string) contains
an acetate ester moiety, i.e. a substructure of the form O–C(=O)–CH3.
"""

from rdkit import Chem

def is_acetate_ester(smiles: str):
    """
    Determines if a molecule given by a SMILES string is (or contains) an acetate ester.
    An acetate ester is defined here as any carboxylic ester where the carboxylic acid component
    is acetic acid (i.e. it contains the acyl group CH3-CO-).
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule contains an acetate ester moiety, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a SMARTS pattern for an acetate ester group.
    # The pattern O-C(=O)[CH3] looks for an oxygen bound to a carbonyl carbon which in turn is bound to a methyl group.
    acetate_smarts = "O-C(=O)[CH3]"
    acetate_pattern = Chem.MolFromSmarts(acetate_smarts)
    if acetate_pattern is None:
        return False, "Failed to create acetate SMARTS pattern"
    
    # Search for the acetate ester substructure in the molecule
    matches = mol.GetSubstructMatches(acetate_pattern)
    if not matches:
        return False, "No acetate ester substructure (O-C(=O)-CH3) found"
    
    return True, "Contains an acetate ester moiety (O-C(=O)-CH3)"

# Example usage (if needed):
if __name__ == "__main__":
    # Test a few examples. For instance: 4-acetoxy benzoic acid has an acetate ester.
    test_smiles = "CC(=O)Oc1ccc(cc1)C(O)=O"
    result, reason = is_acetate_ester(test_smiles)
    print(result, ":", reason)