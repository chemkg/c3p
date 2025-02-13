"""
Classifies: CHEBI:24043 flavones
"""
"""
Classifies: Flavones (compounds with a 2-aryl-1-benzopyran-4-one scaffold and substituted derivatives)
"""
from rdkit import Chem

def is_flavones(smiles: str):
    """
    Determines if a molecule is a flavone (or derivative) based on its SMILES string.
    A flavone is defined by a 2-aryl-1-benzopyran-4-one core.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule contains the flavone core, False otherwise
        str: Reason for classification
    """
    # Parse the input SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a SMARTS pattern for the flavone core.
    # This pattern looks for an aromatic ring attached to a benzopyran-4-one system.
    # The B-ring is represented as "c1ccc(cc1)" and the benzopyran-4-one as "c2oc3ccccc3c(=O)c2".
    flavone_pattern = Chem.MolFromSmarts("c1ccc(cc1)-c2oc3ccccc3c(=O)c2")
    if flavone_pattern is None:
        return False, "Error in SMARTS pattern creation for flavone core"

    # Check if the molecule has a match for the flavone core substructure.
    if not mol.HasSubstructMatch(flavone_pattern):
        return False, "Flavone core (2-aryl-1-benzopyran-4-one) not found in the molecule"
        
    return True, "Molecule contains the flavone core (2-aryl-1-benzopyran-4-one) skeleton"

# Example usage (for testing purposes):
if __name__ == "__main__":
    # Example: diosmetin from user provided examples
    test_smiles = "COC1=CC=C(C=C1O)C1=CC(=O)C2=C(O)C=C(O)C=C2O1"
    result, reason = is_flavones(test_smiles)
    print(f"Result: {result}\nReason: {reason}")