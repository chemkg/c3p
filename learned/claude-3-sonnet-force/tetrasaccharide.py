"""
Classifies: CHEBI:50126 tetrasaccharide
"""
"""
Classifies: CHEBI:36299 tetrasaccharide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, Descriptors
from rdkit import rdBase

def is_tetrasaccharide(smiles: str):
    """
    Determines if a molecule is a tetrasaccharide based on its SMILES string.
    A tetrasaccharide is an oligosaccharide composed of four monosaccharide units.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tetrasaccharide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Count rings and get ring sizes
    ring_info = mol.GetRingInfo()
    num_rings = ring_info.NumRings()
    ring_sizes = [len(ring) for ring in ring_info.AtomRings()]
    
    # Tetrasaccharides should have 4 rings of size 5 or 6
    if num_rings != 4:
        return False, f"Found {num_rings} rings, tetrasaccharides should have 4"
    if any(size < 5 or size > 6 for size in ring_sizes):
        return False, "Found rings of invalid size (should be 5 or 6 membered)"
    
    # Check for glycosidic bonds
    glycosidic_pattern = Chem.MolFromSmarts("[OX2H0][CX4]")
    glycosidic_bonds = mol.GetSubstructMatches(glycosidic_pattern)
    if len(glycosidic_bonds) < 3:
        return False, "Found fewer than 3 glycosidic bonds"
    
    # Count oxygen and carbon atoms
    num_oxygen = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    num_carbon = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    
    # Tetrasaccharides typically have 14-28 carbons and 8-14 oxygens
    if num_carbon < 14 or num_carbon > 28:
        return False, "Carbon count outside typical range for tetrasaccharide"
    if num_oxygen < 8 or num_oxygen > 14:
        return False, "Oxygen count outside typical range for tetrasaccharide"
    
    # Check molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500 or mol_wt > 1000:
        return False, "Molecular weight outside typical range for tetrasaccharide"
    
    return True, "Contains 4 sugar rings linked by glycosidic bonds"