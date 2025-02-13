"""
Classifies: CHEBI:64583 sphingomyelin
"""
"""
Classifies: CHEBI:18023 sphingomyelin
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_sphingomyelin(smiles: str):
    """
    Determines if a molecule is a sphingomyelin based on its SMILES string.
    A sphingomyelin is a phospholipid with a sphingoid base linked to a fatty acid
    via an amide bond, and a phosphorylcholine group ester-linked to the terminal hydroxy group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sphingomyelin, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for sphingoid base pattern
    sphingoid_base_pattern = Chem.MolFromSmarts("[N;X3;H1,H2][C;X4;H1,H2][C;X4;H1,H2][C;X4;H1,H2][C;X4;H1,H2][C;X4;H1,H2][C;X4;H1,H2][C;X4;H1,H2][C;X4;H1,H2][C;X4;H1,H2][C;X4;H1,H2][C;X4;H1,H2][C;X4;H1,H2]")
    if not mol.HasSubstructMatch(sphingoid_base_pattern):
        return False, "No sphingoid base found"

    # Look for fatty acid chain linked via amide bond
    fatty_acid_pattern = Chem.MolFromSmarts("[N;X3;H1,H2][C;X3;H1,H2](=[O;X1])[C;X4;H1,H2]")
    if not mol.HasSubstructMatch(fatty_acid_pattern):
        return False, "No fatty acid chain linked via amide bond"

    # Look for phosphorylcholine group ester-linked to the terminal hydroxy group
    phosphocholine_pattern = Chem.MolFromSmarts("[O;X2][P;X4]([O;X1-])([O;X1-])[O;X2][C;X4][N+;X4]([C;X4])([C;X4])[C;X4]")
    if not mol.HasSubstructMatch(phosphocholine_pattern):
        return False, "No phosphorylcholine group ester-linked to terminal hydroxy group"

    # Check for correct connectivity and positioning
    connectivity_pattern = Chem.MolFromSmarts("[N;X3;H1,H2][C;X4;H1,H2][C;X4;H1,H2][C;X4;H1,H2][C;X4;H1,H2][C;X4;H1,H2][C;X4;H1,H2][C;X4;H1,H2][C;X4;H1,H2][C;X4;H1,H2][C;X4;H1,H2][C;X4;H1,H2][C;X4;H1,H2][C;X4;H1,H2][O;X2][C;X4;H1,H2][N;X3;H1,H2][C;X3;H1,H2](=[O;X1])[C;X4;H1,H2]")
    if not mol.HasSubstructMatch(connectivity_pattern):
        return False, "Key components not correctly connected or positioned"

    # Check for trans double bond in sphingoid base
    trans_double_bond_pattern = Chem.MolFromSmarts("[C;X3]=[C;X3]")
    if not mol.HasSubstructMatch(trans_double_bond_pattern):
        return False, "No trans double bond in sphingoid base"

    # Check for specific stereochemistry of sphingoid base
    stereochemistry_pattern = Chem.MolFromSmarts("[C@H]([N])[C@@H]([O])[C@@H]([C])[C@@H]([C])[C@@H]([C])[C@@H]([C])[C@@H]([C])[C@@H]([C])[C@@H]([C])[C@@H]([C])[C@@H]([C])[C@@H]([C])[O]")
    if not mol.HasSubstructMatch(stereochemistry_pattern):
        return False, "Incorrect stereochemistry of sphingoid base"

    # Check molecular weight and atom counts
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500 or mol_wt > 1000:
        return False, "Molecular weight outside typical range for sphingomyelins"

    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 30 or c_count > 60:
        return False, "Carbon count outside typical range for sphingomyelins"

    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 5 or o_count > 10:
        return False, "Oxygen count outside typical range for sphingomyelins"

    return True, "Contains sphingoid base, fatty acid chain linked via amide bond, and phosphorylcholine group ester-linked to terminal hydroxy group"