"""
Classifies: CHEBI:24654 hydroxy fatty acid
"""
"""
Classifies: CHEBI:36676 hydroxy fatty acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is a hydroxy fatty acid based on its SMILES string.
    A hydroxy fatty acid is any fatty acid carrying one or more hydroxy substituents.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hydroxy fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for carboxylic acid group (-C(=O)O)
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No carboxylic acid group found"

    # Look for hydroxy group (-OH) attached to an aliphatic carbon
    hydroxy_pattern = Chem.MolFromSmarts("[OX2H][CX4]")
    hydroxy_matches = mol.GetSubstructMatches(hydroxy_pattern)
    if len(hydroxy_matches) == 0:
        return False, "No hydroxy group found"

    # Look for aliphatic carbon chain (C~C~C~C~C) connected to the carboxylic acid group
    chain_pattern = Chem.MolFromSmarts("C(=O)O[CX4]~[CX4]~[CX4]~[CX4]~[CX4]")
    chain_matches = mol.GetSubstructMatches(chain_pattern)
    if len(chain_matches) == 0:
        return False, "No aliphatic carbon chain found"

    # Check for fatty acid chain length (at least 4 carbons, excluding the carboxyl carbon)
    n_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6 and atom.GetSymbol() != "C")
    if n_carbons < 4:
        return False, "Carbon chain too short for fatty acid"

    # Exclude molecules with aromatic rings or heterocycles
    if mol.HasSubstructMatch(Chem.MolFromSmarts("a")) or mol.HasSubstructMatch(Chem.MolFromSmarts("*[!#6;!#1;!#8]")):
        return False, "Molecule contains aromatic rings or heterocycles"

    return True, "Molecule contains a carboxylic acid group, at least one hydroxy group, and an aliphatic carbon chain"