"""
Classifies: CHEBI:35757 monocarboxylic acid anion
"""
"""
Classifies: Monocarboxylic acid anion 
Definition: A carboxylic acid anion formed when the carboxy group of a monocarboxylic acid is deprotonated.
"""
from rdkit import Chem

def is_monocarboxylic_acid_anion(smiles: str):
    """
    Determines if a molecule is a monocarboxylic acid anion based on its SMILES string.
    This requires the molecule to contain exactly one deprotonated carboxylic acid group,
    recognized by the substructure [O-]C(=O)[#6].

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is classified as a monocarboxylic acid anion, False otherwise.
        str: Reason for classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the SMARTS pattern for a deprotonated carboxylic acid group (carboxylate)
    carboxylate_pattern = Chem.MolFromSmarts("[O-]C(=O)[#6]")
    if carboxylate_pattern is None:
        return False, "SMARTS pattern creation failure"
    
    # Find all matches of the pattern in the molecule
    matches = mol.GetSubstructMatches(carboxylate_pattern)
    num_matches = len(matches)
    
    # Evaluate if there is exactly one deprotonated carboxylate group
    if num_matches == 0:
        return False, "No deprotonated carboxylate (COO-) group found"
    elif num_matches > 1:
        return False, f"Found {num_matches} carboxylate groups; expected exactly one for a monocarboxylic acid anion"
    
    # Passed all checks: this is a monocarboxylic acid anion
    return True, "Contains exactly one deprotonated carboxylate (COO-) group, consistent with a monocarboxylic acid anion"