"""
Classifies: CHEBI:17088 monoacyl-sn-glycerol 3-phosphate
"""
"""
Classifies: CHEBI:57980 monoacyl-sn-glycerol 3-phosphate
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_monoacyl_sn_glycerol_3_phosphate(smiles: str):
    """
    Determines if a molecule is a monoacyl-sn-glycerol 3-phosphate based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monoacyl-sn-glycerol 3-phosphate, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Create and verify SMARTS patterns
    patterns = {
        'free_phosphate': '[OX2H,O-]P(=O)([OX2H,O-])[OX2]',  # Phosphate with at least 2 free O's
        'ester': '[CX3](=O)[OX2]',
        'glycerol_backbone': '[CH2X4][CHX4][CH2X4]',
        'hydroxyl': '[OX2H]',
        'ethanolamine': 'OCCN',  # Pattern to exclude
        'choline': 'OCC[N+]',    # Pattern to exclude
        'ether': '[CH2,CH][OX2][CH2]C'  # Pattern to exclude ether lipids
    }
    
    compiled_patterns = {}
    for name, pattern in patterns.items():
        pat = Chem.MolFromSmarts(pattern)
        if pat is None:
            return False, f"Invalid SMARTS pattern for {name}"
        compiled_patterns[name] = pat

    # Exclude molecules with ethanolamine or choline groups
    if mol.HasSubstructMatch(compiled_patterns['ethanolamine']):
        return False, "Contains ethanolamine group - not a monoacyl-sn-glycerol 3-phosphate"
    if mol.HasSubstructMatch(compiled_patterns['choline']):
        return False, "Contains choline group - not a monoacyl-sn-glycerol 3-phosphate"

    # Check for free phosphate group
    if not mol.HasSubstructMatch(compiled_patterns['free_phosphate']):
        return False, "No free phosphate group found"

    # Check for single ester group (acyl chain)
    ester_matches = mol.GetSubstructMatches(compiled_patterns['ester'])
    if len(ester_matches) != 1:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 1"

    # Check for glycerol backbone
    if not mol.HasSubstructMatch(compiled_patterns['glycerol_backbone']):
        return False, "No glycerol backbone found"

    # Exclude ether lipids
    if mol.HasSubstructMatch(compiled_patterns['ether']):
        return False, "Contains ether linkage - not a monoacyl-sn-glycerol 3-phosphate"

    # Count oxygen atoms - should be 6 (3 from phosphate, 2 from ester, 1 from hydroxyl)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count != 6:
        return False, f"Incorrect number of oxygen atoms ({o_count}), expected 6"

    # Count phosphorus atoms - should be exactly 1
    p_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
    if p_count != 1:
        return False, f"Found {p_count} phosphorus atoms, need exactly 1"

    # Verify the presence of exactly one hydroxyl group
    hydroxyl_matches = mol.GetSubstructMatches(compiled_patterns['hydroxyl'])
    if len(hydroxyl_matches) != 1:
        return False, f"Found {len(hydroxyl_matches)} free hydroxyl groups, need exactly 1"

    # Additional check for acyl chain length (at least 4 carbons)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 7:  # minimum 3 for glycerol + 4 for shortest acyl chain
        return False, "Acyl chain too short"

    return True, "Contains glycerol backbone with one acyl chain and free phosphate group"