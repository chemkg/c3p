"""
Classifies: CHEBI:24043 flavones
"""
"""
Classifies: Flavones – compounds containing the 2-aryl-1-benzopyran-4-one core.
A flavone is defined as a flavonoid having a 2-phenyl-4H-1-benzopyran-4-one scaffold.
This version uses an improved SMARTS pattern that uses fully aromatic notation and relaxes the match criteria.
"""

from rdkit import Chem

def is_flavones(smiles: str):
    """
    Determines if a molecule is a flavone (or a substituted derivative) based on its SMILES string.
    In this implementation, a molecule is classified as a flavone if it contains a fused 2-aryl-1-benzopyran-4-one core.
    
    The improved SMARTS pattern used is:
      "c1ccc(cc1)-c2oc(=O)c3ccccc23"
    which encodes:
      - a phenyl ring (c1ccc(cc1))
      - attached via a bond (the dash) to a pyranone ring system (c2oc(=O)c3ccccc23)
    This pattern is designed to capture the core structure even when additional substituents are present.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule contains the flavone core, False otherwise.
        str: Explanation for the decision.
    """
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define an improved SMARTS pattern for the flavone core using fully aromatic notation.
    # Note: Using "c" for aromatic atoms throughout ensures more consistent matching.
    flavone_smarts = "c1ccc(cc1)-c2oc(=O)c3ccccc23"
    query = Chem.MolFromSmarts(flavone_smarts)
    if query is None:
        return False, "Error creating SMARTS pattern for flavone core"
    
    # Search for the flavone substructure in the molecule.
    # We do not restrict the match to be of the same size as the query to allow for substituted derivatives.
    if mol.HasSubstructMatch(query):
        return True, "Molecule contains the flavone core (2-aryl-1-benzopyran-4-one) skeleton"
    else:
        return False, "Flavone core (2-aryl-1-benzopyran-4-one) not found in the molecule"

# Example usage for testing:
if __name__ == "__main__":
    # Test with diosmetin, whose SMILES is a known flavone.
    test_smiles = "COC1=CC=C(C=C1O)C1=CC(=O)C2=C(O)C=C(O)C=C2O1"
    result, reason = is_flavones(test_smiles)
    print(f"Result: {result}\nReason: {reason}")