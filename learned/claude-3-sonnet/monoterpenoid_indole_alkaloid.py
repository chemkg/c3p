"""
Classifies: CHEBI:65323 monoterpenoid indole alkaloid
"""
"""
Classifies: CHEBI:51979 monoterpenoid indole alkaloid

A terpenoid indole alkaloid which is biosynthesised from L-tryptophan and 
diisoprenoid (usually secolaganin) building blocks.
"""

from rdkit import Chem

def is_monoterpenoid_indole_alkaloid(smiles: str):
    """
    Determines if a molecule is a monoterpenoid indole alkaloid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monoterpenoid indole alkaloid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for indole substructure
    indole_pattern = Chem.MolFromSmarts("c1ccc2c(c1)nc(c2)C")
    if not mol.HasSubstructMatch(indole_pattern):
        return False, "No indole substructure found"
    
    # Look for terpene substructure (isoprene units)
    terpene_pattern = Chem.MolFromSmarts("[C@](C)(CC=C)")
    if not mol.HasSubstructMatch(terpene_pattern):
        return False, "No terpene substructure found"
    
    # Check for diisoprenoid chain (C10 or C15)
    c10_pattern = Chem.MolFromSmarts("[C@](C)(CC=C)(CC=C)"
    c15_pattern = Chem.MolFromSmarts("[C@](C)(CC=C)(CC=C)(CC=C)")
    if not (mol.HasSubstructMatch(c10_pattern) or mol.HasSubstructMatch(c15_pattern)):
        return False, "No diisoprenoid chain found (C10 or C15)"
    
    # Check for presence of nitrogen
    if mol.GetAtomWithIdx(0).GetAtomicNum() != 7:
        return False, "No nitrogen atom found"
    
    # Check for ring fusion
    rings = mol.GetRingInfo().AtomRings()
    if len(rings) < 2:
        return False, "Not enough fused rings"
    
    return True, "Contains indole nucleus and terpenoid (C10 or C15) backbone"