"""
Classifies: CHEBI:16749 1-phosphatidyl-1D-myo-inositol
"""
from rdkit import Chem

def is_1_phosphatidyl_1D_myo_inositol(smiles: str):
    """
    Determines if a molecule is a 1-phosphatidyl-1D-myo-inositol based on its SMILES string.
    The molecule should have a 1D-myo-inositol with a phosphatidyl group at position 1.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 1-phosphatidyl-1D-myo-inositol, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # SMARTS for 1D-myo-inositol with correct stereochemistry
    inositol_smarts = "[C@H]1(O)[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O"
    inositol_pattern = Chem.MolFromSmarts(inositol_smarts)
    if inositol_pattern is None or not mol.HasSubstructMatch(inositol_pattern):
        return False, "No 1D-myo-inositol moiety found"

    # SMARTS for phosphate group attached at the 1-position
    phosphate_smarts = "O[P](=O)(O)OC1"
    phosphate_linkage_pattern = Chem.MolFromSmarts(phosphate_smarts)
    if phosphate_linkage_pattern is None or not mol.HasSubstructMatch(phosphate_linkage_pattern):
        return False, "No phosphate ester linkage at position 1 of inositol"

    # SMARTS for glycerol backbone linked to the phosphate
    glycerol_smarts = "OCCOP(=O)"
    glycerol_pattern = Chem.MolFromSmarts(glycerol_smarts)
    if glycerol_pattern is None or not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone for phosphatidyl linkage"

    return True, "Contains 1-phosphatidyl-1D-myo-inositol moiety"