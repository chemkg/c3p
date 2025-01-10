"""
Classifies: CHEBI:35366 fatty acid
"""
"""
Classifies: CHEBI:35748 fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_fatty_acid(smiles: str):
    """
    Determines if a molecule is a fatty acid based on its SMILES string.
    A fatty acid is an aliphatic monocarboxylic acid with a chain of 4 to 28 carbons.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group (-COOH)
    carboxylic_acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # Count the number of carbon atoms in the molecule
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 4 or c_count > 28:
        return False, f"Carbon count {c_count} is outside the typical range for fatty acids (4-28)"

    # Check if the molecule is aliphatic (no aromatic rings)
    aromatic_atoms = [atom for atom in mol.GetAtoms() if atom.GetIsAromatic()]
    if aromatic_atoms:
        return False, "Molecule contains aromatic rings, which are not typical for fatty acids"

    # Check for long aliphatic chain (at least 4 carbons in a row)
    aliphatic_chain_pattern = Chem.MolFromSmarts("[CX4][CX4][CX4][CX4]")
    if not mol.HasSubstructMatch(aliphatic_chain_pattern):
        return False, "No aliphatic chain of sufficient length found"

    # Check for unsaturation (double or triple bonds)
    unsaturation_pattern = Chem.MolFromSmarts("[CX3]=[CX3]|[CX2]#[CX2]")
    unsaturation_matches = mol.GetSubstructMatches(unsaturation_pattern)
    if not unsaturation_matches:
        return True, "Saturated fatty acid"
    else:
        return True, "Unsaturated fatty acid"

    return True, "Contains a carboxylic acid group and a long aliphatic chain"