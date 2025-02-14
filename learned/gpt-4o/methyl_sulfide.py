"""
Classifies: CHEBI:86315 methyl sulfide
"""
from rdkit import Chem

def is_methyl_sulfide(smiles: str):
    """
    Determines if a molecule is a methyl sulfide based on its SMILES string.
    A methyl sulfide is defined as a sulfide where at least one of the groups attached to the sulfur is a methyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a methyl sulfide, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the pattern for methyl sulfide (S-methyl)
    methyl_sulfide_pattern = Chem.MolFromSmarts("CS")

    # Check if the molecule has a sulfur atom with a bonded methyl group
    if mol.HasSubstructMatch(methyl_sulfide_pattern):
        return True, "Contains sulfur atom directly bonded to a methyl group"

    return False, "Does not have a sulfur atom bonded to a methyl group"

# Example usages of the function
example_smiles = [
    "O=C(O)[C@@H](N)CCCCSC",
    "CSc1nc(N)nc(N)n1",
    "O=C(O)[C@@H](N)CCCCCSC",
    "C(=C/C(O)=O)(\CCCSC)/C(=O)O",
    "CSCCCO",
    "CCSC"
]

for smiles in example_smiles:
    result, reason = is_methyl_sulfide(smiles)
    print(f"SMILES: {smiles} -> Is Methyl Sulfide: {result}, Reason: {reason}")