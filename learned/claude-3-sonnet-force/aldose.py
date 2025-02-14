"""
Classifies: CHEBI:15693 aldose
"""
"""
Classifies: CHEBI:17395 aldose
An aldose is an aldehydic parent sugar (polyhydroxy aldehyde H[CH(OH)]nC(=O)H, n >= 2) 
or its intramolecular hemiacetal.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_aldose(smiles: str):
    """
    Determines if a molecule is an aldose based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aldose, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Count aldehyde groups (C=O)
    aldehyde_pattern = Chem.MolFromSmarts("C=O")
    aldehyde_matches = mol.GetSubstructMatches(aldehyde_pattern)
    n_aldehydes = len(aldehyde_matches)
    
    # Count hydroxyl groups (-OH)
    hydroxyl_pattern = Chem.MolFromSmarts("[OH]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    n_hydroxyls = len(hydroxyl_matches)
    
    # Count carbon atoms
    n_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    
    # Aldoses must have at least 1 aldehyde, 2 hydroxyls, and a minimum number of carbon atoms
    if n_aldehydes == 0 or n_hydroxyls < 2 or n_carbons < 4:
        return False, "Not enough aldehyde/hydroxyl groups or carbon atoms for aldose"
    
    # Check for linear aldose (H[CH(OH)]nC(=O)H or C(=O)H[CH(OH)]n)
    linear_aldose_pattern = Chem.MolFromSmarts("[CH](O)[CH](O)[CH](O)[CH](=O)O")
    if mol.HasSubstructMatch(linear_aldose_pattern):
        return True, "Linear aldose detected"
    
    # Check for cyclic hemiacetal
    ring_info = mol.GetRingInfo()
    for ring in ring_info.AtomRings():
        # Check if ring contains aldehyde and hydroxyl groups
        ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
        ring_aldehydes = sum(1 for atom in ring_atoms if atom.GetAtomicNum() == 8 and atom.GetDegree() == 1)
        ring_hydroxyls = sum(1 for atom in ring_atoms if atom.GetAtomicNum() == 8 and atom.GetDegree() == 2)
        if ring_aldehydes == 1 and ring_hydroxyls >= 1:
            return True, "Cyclic hemiacetal aldose detected"
    
    # Check for molecular weight (typical aldoses are >120 Da)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 120:
        return False, "Molecular weight too low for aldose"
    
    # Check for carboxylic acid groups (-COOH) or other functional groups
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)O")
    if mol.HasSubstructMatch(carboxyl_pattern):
        return False, "Contains carboxylic acid group, not an aldose"
    
    return False, "No aldose patterns found"