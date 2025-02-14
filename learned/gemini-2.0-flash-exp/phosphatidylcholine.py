"""
Classifies: CHEBI:64482 phosphatidylcholine
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_phosphatidylcholine(smiles: str):
    """
    Determines if a molecule is a phosphatidylcholine based on its SMILES string.
    A phosphatidylcholine is a glycerol backbone with a phosphate group attached to the
    third carbon, a choline head group and two fatty acid chains attached via ester bonds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phosphatidylcholine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Glycerol Backbone (correct stereochemistry)
    glycerol_pattern = Chem.MolFromSmarts("[C@H]([OX2])[C]([OX2])[CH2][OX2]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"
    
    # 2. Phosphate Group and Choline Head Group
    phosphate_choline_pattern = Chem.MolFromSmarts("[CH2][O][P](=[O])([O-])OCC[N+](C)(C)C")
    if not mol.HasSubstructMatch(phosphate_choline_pattern):
        return False, "No phosphate and/or choline head group found"

    # 3. Two ester groups attached to the glycerol at positions 1 and 2 and the complete structure
    complete_pattern = Chem.MolFromSmarts("[C@H]([OX2][CX3](=[OX1]))([OX2][CX3](=[OX1]))[CH2X4][OX2][PX4](=[OX1])([OX2-])OCC[N+](C)(C)C")
    if not mol.HasSubstructMatch(complete_pattern):
         return False, "Ester groups are not on the 1 and 2 positions of the glycerol backbone and/or missing components"

    # Check for fatty acid chains (long carbon chains attached to esters)
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if len(fatty_acid_matches) < 2:
        return False, f"Missing fatty acid chains, got {len(fatty_acid_matches)}"

    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    p_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)

    if c_count < 8:
        return False, "Too few carbons for a phosphatidylcholine"

    if o_count < 6:
        return False, "Must have at least 6 oxygens"
    
    if n_count != 1:
        return False, "Must have exactly one nitrogen (choline)"
    
    if p_count != 1:
        return False, "Must have exactly one phosphorus (phosphate)"
    
    return True, "Meets all criteria for phosphatidylcholine"