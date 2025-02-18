"""
Classifies: CHEBI:28494 cardiolipin
"""
"""
Classifies: CHEBI:28425 cardiolipin
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_cardiolipin(smiles: str):
    """
    Determines if a molecule is a cardiolipin based on its SMILES string.
    A cardiolipin is a phosphatidylglycerol composed of two phosphatidic acid molecules linked to a central glycerol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cardiolipin, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for exactly two phosphorus atoms (from phosphate groups)
    p_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
    if p_count != 2:
        return False, f"Found {p_count} phosphorus atoms, need exactly 2"

    # Detect ester groups (R-O-C=O)
    ester_pattern = Chem.MolFromSmarts("[OX2]C(=O)[#6]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 4:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 4"

    # Detect glycerol backbones (three carbons with oxygen substituents)
    # Updated pattern: Each carbon must have at least one oxygen (any type)
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4]")
    glycerol_matches = mol.GetSubstructMatches(glycerol_pattern)
    
    # Filter matches to require oxygen on each carbon
    valid_glycerols = []
    for match in glycerol_matches:
        c1, c2, c3 = match
        if (mol.GetAtomWithIdx(c1).GetDegree() == 4 and
            mol.GetAtomWithIdx(c2).GetDegree() == 4 and
            mol.GetAtomWithIdx(c3).GetDegree() == 4 and
            any(a.GetAtomicNum() == 8 for a in mol.GetAtomWithIdx(c1).GetNeighbors()) and
            any(a.GetAtomicNum() == 8 for a in mol.GetAtomWithIdx(c2).GetNeighbors()) and
            any(a.GetAtomicNum() == 8 for a in mol.GetAtomWithIdx(c3).GetNeighbors())):
            valid_glycerols.append(match)
    
    if len(valid_glycerols) < 3:
        return False, f"Found {len(valid_glycerols)} valid glycerol backbones, need at least 3"

    # Check phosphate groups (O-P-O linkages)
    phosphate_pattern = Chem.MolFromSmarts("[OX2]P(=O)([OX2])[OX2]")
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if len(phosphate_matches) < 2:
        return False, f"Found {len(phosphate_matches)} phosphate groups, need at least 2"

    # Verify molecular weight (typical cardiolipins >1200 Da)
    mol_wt = Descriptors.ExactMolWt(mol)
    if mol_wt < 1200:
        return False, f"Molecular weight too low ({mol_wt:.1f} Da)"

    # Check fatty acid chain lengths (at least 12 carbons)
    chain_pattern = Chem.MolFromSmarts("[#6]C(=O)O[#6]")
    for match in mol.GetSubstructMatches(chain_pattern):
        chain_start = match[0]
        atom = mol.GetAtomWithIdx(chain_start)
        chain_length = 0
        
        # Traverse the chain away from carbonyl
        while True:
            neighbors = [n for n in atom.GetNeighbors() 
                        if n.GetAtomicNum() == 6 
                        and n.GetIdx() != match[1]]  # Avoid carbonyl
            if len(neighbors) != 1:
                break
            chain_length += 1
            atom = neighbors[0]
            
        if chain_length < 12:
            return False, f"Fatty acid chain too short ({chain_length} carbons)"

    return True, "Contains two phosphatidic acids linked via central glycerol with four ester groups and two phosphates"