"""
Classifies: CHEBI:27933 beta-lactam antibiotic
"""
"""
Classifies: CHEBI:35627 beta-lactam antibiotic
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_beta_lactam_antibiotic(smiles: str):
    """
    Determines if a molecule is a beta-lactam antibiotic based on its SMILES string.
    A beta-lactam antibiotic must contain a beta-lactam ring (a four-membered cyclic amide)
    and typically has antibiotic properties.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a beta-lactam antibiotic, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a more specific beta-lactam ring pattern (four-membered cyclic amide)
    beta_lactam_pattern = Chem.MolFromSmarts("[C,c]1(=O)[N,n][C,c][C,c]1")
    if not mol.HasSubstructMatch(beta_lactam_pattern):
        return False, "No beta-lactam ring found"

    # Check for antibiotic-like properties (e.g., presence of nitrogen and oxygen atoms)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if n_count < 1:
        return False, "No nitrogen atoms found, unlikely to be an antibiotic"
    if o_count < 1:
        return False, "No oxygen atoms found, unlikely to be an antibiotic"

    # Check molecular weight (beta-lactam antibiotics typically have MW > 200 Da)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200:
        return False, "Molecular weight too low for a beta-lactam antibiotic"

    # Check for common antibiotic functional groups (e.g., amide, carboxyl, etc.)
    # Relaxed criteria to avoid false negatives
    amide_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[NX3H2]")
    carboxyl_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2H1]")
    ester_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2][CX4]")
    
    if not (mol.HasSubstructMatch(amide_pattern) or 
            mol.HasSubstructMatch(carboxyl_pattern) or 
            mol.HasSubstructMatch(ester_pattern)):
        # Allow molecules without these groups if they have other antibiotic-like features
        # Check for other common antibiotic substructures
        antibiotic_substructure = Chem.MolFromSmarts("[C,c]1(=O)[N,n][C,c][C,c]1")
        if not mol.HasSubstructMatch(antibiotic_substructure):
            return False, "No common antibiotic functional groups or substructures found"

    return True, "Contains a beta-lactam ring and exhibits antibiotic-like properties"