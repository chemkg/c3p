"""
Classifies: CHEBI:16337 phosphatidic acid
"""
from rdkit import Chem

def is_phosphatidic_acid(smiles: str):
    """
    Determines if a molecule is a phosphatidic acid based on its SMILES string.
    Phosphatidic acids have a glycerol backbone with two fatty acid chains and a phosphate group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phosphatidic acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Improved glycerol backbone pattern including flexible chirality configurations
    glycerol_pattern = Chem.MolFromSmarts("O[C@@H](CO)C(O)CO")  # Considering flexible SRC patterns
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone with required connectivity found"

    # Look for exactly 2 ester linkages including the carbonyl group
    ester_pattern = Chem.MolFromSmarts("C(=O)O[CH]")  # Flexible anchoring to oxygen and carbon
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 2:
        return False, f"Found {len(ester_matches)} ester linkages, need exactly 2"
    
    # Confirm phosphate group presence with correct oxygen attachments
    phosphate_pattern = Chem.MolFromSmarts("P(=O)(O)O")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate group with necessary attachments found"
    
    return True, "Contains glycerol backbone with 2 fatty acid chains and a phosphate group"

# Example usage
smiles = "P(OC[C@H](OC(=O)CCCCCCCCCCCCCCCC)COC(=O)CCCCCCCCCCCCC)(O)(O)=O"
result, reason = is_phosphatidic_acid(smiles)
print(f"Is phosphatidic acid: {result}, Reason: {reason}")