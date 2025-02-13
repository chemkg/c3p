"""
Classifies: CHEBI:27093 tricarboxylic acid
"""
"""
Classifies: Tricarboxylic acid – An oxoacid containing three carboxy groups.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_tricarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a tricarboxylic acid based on its SMILES string.
    
    A tricarboxylic acid is defined as an oxoacid containing three carboxy groups.
    Carboxyl groups are identified here by the substructure C(=O)O in the protonated form
    or C(=O)[O-] in the deprotonated form.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule contains exactly three carboxyl groups, False otherwise
        str: Reason for the classification
    """
    # Parse the SMILES string into an RDKit Mol object.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS patterns for the carboxyl group in both protonated and deprotonated forms.
    # Pattern for carboxylic acid: C(=O)[OH]
    carboxyl_neutral = Chem.MolFromSmarts("C(=O)[OX2H1]")
    # Pattern for carboxylate anion: C(=O)[O-]
    carboxyl_anion = Chem.MolFromSmarts("C(=O)[OX1-]")
    
    # Get matches for each pattern.
    matches_neutral = set(mol.GetSubstructMatches(carboxyl_neutral))
    matches_anion = set(mol.GetSubstructMatches(carboxyl_anion))
    
    # Create a union of unique matches.
    carboxyl_matches = matches_neutral.union(matches_anion)
    n_carboxyl = len(carboxyl_matches)
    
    # Check for the presence of exactly three carboxyl groups.
    if n_carboxyl != 3:
        return False, f"Found {n_carboxyl} carboxyl group(s), but exactly 3 are required."
    
    return True, "Molecule contains exactly 3 carboxyl groups and is classified as a tricarboxylic acid."

# Example usage (uncomment for testing):
# test_smiles = "OC(=O)CC(O)(CC(O)=O)C(O)=O"  # Example: citric acid has three carboxyl groups.
# result, reason = is_tricarboxylic_acid(test_smiles)
# print(result, reason)