"""
Classifies: CHEBI:65323 monoterpenoid indole alkaloid
"""
"""
Classifies: monoterpenoid indole alkaloids
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_monoterpenoid_indole_alkaloid(smiles: str):
    """
    Determines if a molecule is a monoterpenoid indole alkaloid based on its SMILES string.
    These compounds are derived from tryptophan and typically secologanin.

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

    # Check molecular weight (expanded range for MIAs)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 250 or mol_wt > 900:
        return False, f"Molecular weight {mol_wt:.1f} outside typical MIA range"

    # Count basic statistics
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    ring_count = rdMolDescriptors.CalcNumRings(mol)

    # Basic requirements
    if n_count == 0:
        return False, "No nitrogen atoms found - required for alkaloid"
    if c_count < 16:
        return False, f"Too few carbons ({c_count}) for MIA structure"
    if ring_count < 4:
        return False, f"Too few rings ({ring_count}) for MIA structure"

    # Look for indole or modified indole cores using multiple patterns
    indole_patterns = [
        "c1ccc2[nH]ccc2c1",  # Basic indole
        "c1ccc2nccc2c1",     # Modified indole
        "c1ccc2N=CCc2c1",    # Another modified form
        "c1ccc2NCCc2c1",     # Dihydroindole
        "C1=CC=C2C(=C1)NC=C2",  # Alternative representation
        "C1=CC=C2C(=C1)N=CC2"   # Another variant
    ]
    
    has_indole = False
    for pattern in indole_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            has_indole = True
            break
            
    if not has_indole:
        return False, "No indole or modified indole core found"

    # Check for nitrogen-containing rings (essential for MIA skeleton)
    n_ring_pattern = Chem.MolFromSmarts("[N;R]")
    if len(mol.GetSubstructMatches(n_ring_pattern)) < 2:
        return False, "Insufficient nitrogen-containing rings"

    # Look for typical MIA structural features
    features = []
    
    # Check for ester groups (common in MIAs)
    if mol.HasSubstructMatch(Chem.MolFromSmarts("C(=O)OC")):
        features.append("ester group")
    
    # Check for ethyl/vinyl groups (common in monoterpene portion)
    if mol.HasSubstructMatch(Chem.MolFromSmarts("CC")) or mol.HasSubstructMatch(Chem.MolFromSmarts("C=C")):
        features.append("ethyl/vinyl group")

    # Count sp2 carbons (common in terpene portion)
    sp2_carbons = len(mol.GetSubstructMatches(Chem.MolFromSmarts("[C]=[C,N,O]")))
    if sp2_carbons > 0:
        features.append(f"{sp2_carbons} sp2 carbons")

    # Look for specific MIA structural patterns
    mia_patterns = [
        "[C]1[C]2[N][C]([C]1)[C]([C]2)",  # Common MIA ring fusion pattern
        "[C]1[C]2[N][C]([C]1)[C]([C]2)[O,N]", # Oxygenated/N-containing variant
        "[C]1[C]2[N][C]([C]1)[C]([C]2)C(=O)", # Carbonyl-containing variant
    ]
    
    mia_pattern_matches = 0
    for pattern in mia_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            mia_pattern_matches += 1

    # Calculate revised complexity score
    complexity_contributors = [
        ring_count * 2,      # Rings
        n_count * 3,         # Nitrogens
        sp2_carbons,         # sp2 carbons
        len(features) * 2,   # Typical features
        mia_pattern_matches * 4  # Specific MIA patterns
    ]
    complexity_score = sum(complexity_contributors)
    
    if complexity_score < 15:
        return False, f"Complexity score ({complexity_score}) too low for typical MIA"

    # Additional structural requirements
    if ring_count > 12 and mia_pattern_matches == 0:
        return False, "Complex structure lacks typical MIA ring patterns"
        
    if n_count > 5 and mia_pattern_matches == 0:
        return False, "High nitrogen count but lacks typical MIA patterns"

    # Check C/N ratio (typically between 8-20 for MIAs)
    cn_ratio = c_count / n_count
    if cn_ratio < 8 or cn_ratio > 20:
        return False, f"C/N ratio ({cn_ratio:.1f}) outside typical MIA range"

    feature_str = ", ".join(features) if features else "typical structural features"
    return True, f"Contains indole core, complex ring system ({ring_count} rings), and {feature_str}"