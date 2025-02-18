"""
Classifies: CHEBI:24026 fatty alcohol
"""
"""
Classifies: CHEBI:30843 fatty alcohol
"""
from rdkit import Chem
from rdkit.Chem import Descriptors

def is_fatty_alcohol(smiles: str):
    """
    Determines if a molecule is a fatty alcohol based on its SMILES string.
    A fatty alcohol is an aliphatic alcohol with a carbon chain length of 3 to 27 atoms.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty alcohol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for at least one alcohol group (-OH)
    alcohol_pattern = Chem.MolFromSmarts("[OX2H]")
    if not mol.HasSubstructMatch(alcohol_pattern):
        return False, "No alcohol group found"

    # Count the number of carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 3 or c_count > 27:
        return False, f"Carbon chain length {c_count} is outside the range of 3 to 27"

    # Ensure the molecule is aliphatic (no aromatic rings)
    aromatic_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetIsAromatic())
    if aromatic_atoms > 0:
        return False, "Molecule contains aromatic rings, not aliphatic"

    # Count the number of oxygen atoms (should be at least one, but can be more)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 1:
        return False, "No oxygen atoms found"

    # Check for unsaturation (double or triple bonds)
    unsaturation = sum(1 for bond in mol.GetBonds() if bond.GetBondType() in [Chem.BondType.DOUBLE, Chem.BondType.TRIPLE])
    if unsaturation > 0:
        # Allow unsaturated fatty alcohols
        pass

    # Check for other functional groups that are not typical in fatty alcohols
    # For example, carboxylic acids, esters, etc.
    other_functional_groups = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")  # Carboxylic acid
    if mol.HasSubstructMatch(other_functional_groups):
        return False, "Molecule contains carboxylic acid groups, not a fatty alcohol"

    # Check if the molecule is primarily a hydrocarbon chain with alcohol groups
    # This is a heuristic and may need refinement
    h_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 1)
    if h_count < 2 * c_count - 2:  # Rough check for hydrocarbon-like structure
        return False, "Molecule does not have a hydrocarbon-like structure"

    return True, f"Aliphatic alcohol with {c_count} carbon atoms and {o_count} alcohol groups"