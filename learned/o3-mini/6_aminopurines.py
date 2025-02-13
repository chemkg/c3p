"""
Classifies: CHEBI:20706 6-aminopurines
"""
"""
Classifies: 6-aminopurines (any compound having 6-aminopurine (adenine) as part of its structure)
"""
from rdkit import Chem

def is_6_aminopurines(smiles: str):
    """
    Determines if a molecule contains a 6-aminopurine (adenine) moiety.
    This function uses a flexible SMARTS pattern representing the adenine core:
    a purine ring system bearing an exocyclic amino group.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule contains a 6-aminopurine (adenine) substructure, False otherwise.
        str: Detailed reason for the classification.
    """
    # Convert the SMILES string to an RDKit molecule object.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a more flexible SMARTS pattern for the adenine (6-aminopurine) core.
    # The pattern "n1c2nc(nc(n2)c1)N" describes an aromatic bicyclic ring (purine)
    # with an exocyclic amino group that corresponds to the 6-aminopurine functionality.
    adenine_smarts = "n1c2nc(nc(n2)c1)N"
    adenine_query = Chem.MolFromSmarts(adenine_smarts)
    if adenine_query is None:
        return False, "Error creating adenine substructure query"
    
    # Check if the molecule matches the adenine SMARTS pattern.
    # We are not enforcing chirality here, so this should tolerate modifications.
    if mol.HasSubstructMatch(adenine_query):
        return True, "Molecule contains the 6-aminopurine (adenine) moiety"
    else:
        return False, "6-aminopurine (adenine) moiety not found in the molecule"

# Example usage (to test the functionality, uncomment below lines):
# test_smiles = "Cn1cnc(N)c2ncnc12"  # This is 3-methyladenine; expected to be True.
# result, reason = is_6_aminopurines(test_smiles)
# print(result, reason)