"""
Classifies: CHEBI:62499 methyl-branched fatty acid
"""
"""
Classifies: CHEBI:27718 Methyl-branched fatty acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_methyl_branched_fatty_acid(smiles: str):
    """
    Determines if a molecule is a methyl-branched fatty acid based on its SMILES string.
    A methyl-branched fatty acid is a branched-chain fatty acid containing methyl branches only.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a methyl-branched fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check if molecule contains only C, H, and O atoms
    allowed_atoms = [6, 1, 8]  # C, H, O
    if any(atom.GetAtomicNum() not in allowed_atoms for atom in mol.GetAtoms()):
        return False, "Molecule contains atoms other than C, H, and O"

    # Check for carboxylic acid group
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No carboxylic acid group found"

    # Check for methyl branches
    methyl_branch_pattern = Chem.MolFromSmarts("CC")
    methyl_branch_matches = mol.GetSubstructMatches(methyl_branch_pattern)
    if not methyl_branch_matches:
        return False, "No methyl branches found"

    # Check if all branching groups are methyl groups
    for atom in mol.GetAtoms():
        if atom.GetDegree() > 2 and atom.GetAtomicNum() == 6:
            neighbors = [mol.GetAtomWithIdx(n.GetIdx()) for n in atom.GetNeighbors()]
            methyl_groups = sum(1 for n in neighbors if n.GetAtomicNum() == 6 and n.GetDegreeWithPattern("CC") == 1)
            if atom.GetDegree() - methyl_groups > 2:
                return False, "Found non-methyl branching groups"

    return True, "Contains a carboxylic acid group and only methyl branches"