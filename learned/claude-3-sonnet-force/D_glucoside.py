"""
Classifies: CHEBI:35436 D-glucoside
"""
"""
Classifies: CHEBI:27815 D-glucoside
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_D_glucoside(smiles: str):
    """
    Determines if a molecule is a D-glucoside based on its SMILES string.
    A D-glucoside is a glucoside in which the glycosidic group is derived from D-glucose.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a D-glucoside, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for D-glucose substructure
    d_glucose_pattern = Chem.MolFromSmarts("OC[C@H]1O[C@@H]([C@H]([C@@H]([C@H]1O)O)O)CO")
    if not mol.HasSubstructMatch(d_glucose_pattern):
        return False, "No D-glucose substructure found"
    
    # Look for pyranose ring
    pyranose_pattern = Chem.MolFromSmarts("O[C@H]1[C@H]([C@H]([C@@H]([C@H]([C@@H]1O)O)O)O)O")
    if not mol.HasSubstructMatch(pyranose_pattern):
        return False, "No pyranose ring found"
    
    # Look for glycosidic bonds (O-C-O)
    glycosidic_bond_pattern = Chem.MolFromSmarts("[OX2][CR][OX2]")
    glycosidic_bond_matches = mol.GetSubstructMatches(glycosidic_bond_pattern)
    if not glycosidic_bond_matches:
        return False, "No glycosidic bonds found"
    
    # Check molecular weight - glucosides typically >300 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:
        return False, "Molecular weight too low for glucoside"

    return True, "Contains D-glucose substructure, pyranose ring, and glycosidic bonds"