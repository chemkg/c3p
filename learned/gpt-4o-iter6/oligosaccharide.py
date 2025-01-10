"""
Classifies: CHEBI:50699 oligosaccharide
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_oligosaccharide(smiles: str):
    """
    Determines if a molecule is an oligosaccharide based on its SMILES string.
    An oligosaccharide is a compound with monosaccharide units joined by glycosidic linkages.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an oligosaccharide, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string into a molecule object
    mol = Chem.MolFromSmiles(smiles)
    
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS patterns for various sugar units (hexoses and pentoses with varying stereochemistry)
    sugar_patterns = [
        Chem.MolFromSmarts("[C@H]1(O)C(O)C(O)C(O)C(O)C1"),  # Just one example hexose pattern
        Chem.MolFromSmarts("[C@@H]1(O)C(O)C(O)C(O)C(O)C1"),  # Stereoisomer of the hexose
        Chem.MolFromSmarts("[C@H]1(O)C(O)C(O)C(O)C1"),  # Example pentose pattern
        Chem.MolFromSmarts("[C@@H]1(O)C(O)C(O)C(O)C1")   # Another pentose pattern
        # More patterns can be added as needed
    ]
    
    # Scan the molecule for matching sugar substructures
    sugar_matches = 0
    for pattern in sugar_patterns:
        sugar_matches += len(mol.GetSubstructMatches(pattern))
    
    if sugar_matches < 2:
        return False, f"Found {sugar_matches} sugar units, need at least 2 for oligosaccharide"
    
    # Define more specific glycosidic linkage pattern (anomeric carbon to hydroxyl)
    glycosidic_linkage_pattern = Chem.MolFromSmarts("[C@H]1(OCCC1)-O-*-C")
    
    # Search for glycosidic linkages
    glycosidic_matches = mol.GetSubstructMatches(glycosidic_linkage_pattern)
    
    if len(glycosidic_matches) < 1:
        return False, "No glycosidic linkages found, important for oligosaccharide classification"
    
    return True, "Contains sufficient monosaccharide units linked by glycosidic bonds"