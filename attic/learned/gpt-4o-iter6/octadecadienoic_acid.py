"""
Classifies: CHEBI:25627 octadecadienoic acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_octadecadienoic_acid(smiles: str):
    """
    Determines if a molecule is an octadecadienoic acid based on its SMILES string.
    An octadecadienoic acid is a straight-chain, C18 polyunsaturated fatty acid with two C=C double bonds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an octadecadienoic acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check total carbon count in the molecule
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    # Allow +/- 1 for variations like functional groups directly on the chain
    if not (17 <= carbon_count <= 19):
        return False, f"Expected around 18 carbons, found {carbon_count}"

    # Check for two C=C double bonds
    double_bond_count = sum(1 for bond in mol.GetBonds() if bond.GetBondType() == Chem.BondType.DOUBLE and
                            bond.GetBeginAtom().GetAtomicNum() == 6 and bond.GetEndAtom().GetAtomicNum() == 6)
    if double_bond_count != 2:
        return False, f"Expected 2 double bonds (C=C), found {double_bond_count}"

    # Check for terminal carboxylic acid group or its anions (-COOH or -COO-)
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    carboxylate_pattern = Chem.MolFromSmarts("C(=O)[O-]")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern) and not mol.HasSubstructMatch(carboxylate_pattern):
        return False, "No terminal carboxylic acid group found"

    return True, "Molecule is a C18 polyunsaturated fatty acid with two C=C double bonds and a terminal carboxylic acid group"