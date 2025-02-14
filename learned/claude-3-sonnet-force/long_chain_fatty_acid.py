"""
Classifies: CHEBI:15904 long-chain fatty acid
"""
"""
Classifies: CHEBI:38804 long-chain fatty acid
A fatty acid with a chain length ranging from C13 to C22.
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_long_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a long-chain fatty acid (C13 to C22) based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a long-chain fatty acid, False otherwise
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
    
    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if o_count != 2:
        return False, "Must have exactly 2 oxygens (carboxyl group)"
    
    # Check molecular formula
    mol_formula = rdMolDescriptors.CalcMolFormula(mol)
    if not mol_formula.startswith("C"):
        return False, "Molecular formula does not start with C"
    
    # Check for aromatic atoms
    if any(atom.GetIsAromatic() for atom in mol.GetAtoms()):
        return False, "Aromatic atoms not allowed"
    
    double_bond_pattern = Chem.MolFromSmarts("=")
    n_double_bonds = len(mol.GetSubstructMatches(double_bond_pattern))
    if n_double_bonds > 4:
        return False, "Too many double bonds (>4)"
    
    # Check carbon chain length (C13 to C22)
    if c_count < 13 or c_count > 22:
        return False, f"Carbon chain length not in range C13 to C22 (found C{c_count})"
    
    return True, "Contains carboxylic acid group and long carbon chain (C13 to C22)"