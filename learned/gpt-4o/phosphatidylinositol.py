"""
Classifies: CHEBI:28874 phosphatidylinositol
"""
from rdkit import Chem

def is_phosphatidylinositol(smiles: str):
    """
    Determines if a molecule is a phosphatidylinositol based on its SMILES string.
    A phosphatidylinositol involves a phosphatidyl group esterified to one of the
    hydroxy groups of inositol, a hexahydroxy cyclohexane.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phosphatidylinositol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES string to RDKit molecule object.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define substructure patterns.
    inositol_pattern = Chem.MolFromSmarts("C1([OH])C([OH])C([OH])C([OH])C([OH])C1[OH]")
    phosphatidyl_pattern = Chem.MolFromSmarts("[O-]P(=O)(OC)OC(C(=O))OC(C(=O))")
    
    # Check for inositol ring presence.
    if not mol.HasSubstructMatch(inositol_pattern):
        return False, "No inositol ring with sufficient hydroxylations found"
    
    # Check for phosphatidyl group presence.
    if not mol.HasSubstructMatch(phosphatidyl_pattern):
        return False, "No phosphatidyl group associated with glycerol found"
    
    # Check for connection between inositol and phosphatidyl group.
    connected_pattern = Chem.MolFromSmarts("C1([OH])C([OH])C([OH])C([OH])C([OH])C1[OH]O[P](=O)(=/[O])/[O][CH2]C(=O)")
    if not mol.HasSubstructMatch(connected_pattern):
        return False, "Inositol and phosphatidyl groups not adequately connected"
    
    return True, "The molecule is a phosphatidylinositol with correctly established linkages"