"""
Classifies: CHEBI:23899 icosanoid
"""
"""
Classifies: CHEBI:36080 icosanoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_icosanoid(smiles: str):
    """
    Determines if a molecule is an icosanoid based on its SMILES string.
    An icosanoid is a signaling molecule derived from the oxidation of C20 essential fatty acids.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an icosanoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for core icosanoid structure (16-22 carbons)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 16 or c_count > 22:
        return False, f"Expected 16-22 carbons, found {c_count}"

    # Check for multiple double bonds (typically 1-8)
    double_bonds = sum(1 for bond in mol.GetBonds() if bond.GetBondType() == Chem.BondType.DOUBLE)
    if double_bonds < 1:
        return False, f"Expected at least 1 double bond, found {double_bonds}"

    # Check for oxygen-containing functional groups (at least 2)
    oxygen_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if oxygen_atoms < 2:
        return False, f"Expected at least 2 oxygen atoms, found {oxygen_atoms}"

    # More specific patterns for characteristic icosanoid structures
    # 1. Prostaglandin-like structure (cyclopentane ring with oxygen-containing groups)
    prostaglandin_pattern = Chem.MolFromSmarts("[C]1[C][C][C][C]1([C]=O)([OH])")
    # 2. Leukotriene-like structure (conjugated triene system)
    leukotriene_pattern = Chem.MolFromSmarts("[C]=[C][C]=[C][C]=[C][C](=O)[OH]")
    # 3. Epoxide structure (common in EETs)
    epoxide_pattern = Chem.MolFromSmarts("[OX2]1[CX4][CX4]1")
    # 4. Hydroperoxide structure (common in HPETEs)
    hydroperoxide_pattern = Chem.MolFromSmarts("[OH][OX2]")

    has_prostaglandin = mol.HasSubstructMatch(prostaglandin_pattern)
    has_leukotriene = mol.HasSubstructMatch(leukotriene_pattern)
    has_epoxide = mol.HasSubstructMatch(epoxide_pattern)
    has_hydroperoxide = mol.HasSubstructMatch(hydroperoxide_pattern)

    if not (has_prostaglandin or has_leukotriene or has_epoxide or has_hydroperoxide):
        return False, "No characteristic icosanoid structural patterns found"

    # Check for specific functional groups (hydroxyl, carboxyl, epoxide, peroxide, ester)
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    carboxyl_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    peroxide_pattern = Chem.MolFromSmarts("[OX2][OX2]")
    ester_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H0]")

    has_hydroxyl = mol.HasSubstructMatch(hydroxyl_pattern)
    has_carboxyl = mol.HasSubstructMatch(carboxyl_pattern)
    has_peroxide = mol.HasSubstructMatch(peroxide_pattern)
    has_ester = mol.HasSubstructMatch(ester_pattern)

    if not (has_hydroxyl or has_carboxyl or has_peroxide or has_ester):
        return False, "No hydroxyl, carboxyl, peroxide, or ester groups found"

    return True, "Contains characteristic icosanoid structure with multiple double bonds and oxygen-containing functional groups"