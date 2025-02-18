"""
Classifies: CHEBI:35681 secondary alcohol
"""
"""
Classifies: Secondary Alcohols (CHEBI:35678)
"""
from rdkit import Chem
from rdkit.Chem import HybridizationType

def is_secondary_alcohol(smiles: str) -> tuple[bool, str]:
    """
    Determines if a molecule is a secondary alcohol based on its SMILES string.
    A secondary alcohol has a hydroxyl group (-OH) attached to a saturated carbon atom
    that is bonded to exactly two other carbon atoms.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a secondary alcohol, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Find all hydroxyl oxygen atoms (-OH groups)
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    
    if not hydroxyl_matches:
        return False, "No hydroxyl group found"

    for match in hydroxyl_matches:
        oxygen_idx = match[0]
        oxygen = mol.GetAtomWithIdx(oxygen_idx)
        
        # Get the connected carbon (should be only one neighbor for -OH)
        neighbors = oxygen.GetNeighbors()
        if len(neighbors) != 1:
            continue  # Not a typical hydroxyl group
        
        carbon = neighbors[0]
        if carbon.GetAtomicNum() != 6:
            continue  # Hydroxyl not attached to carbon
            
        # Check for sp³ hybridization (saturated carbon)
        if carbon.GetHybridization() != HybridizationType.SP3:
            continue  # Exclude aromatic/unsaturated carbons
            
        # Count carbon neighbors (must have exactly two)
        carbon_neighbors = [nbr for nbr in carbon.GetNeighbors() 
                          if nbr.GetAtomicNum() == 6 
                          and nbr.GetIdx() != oxygen_idx]
        
        if len(carbon_neighbors) == 2:
            return True, "Hydroxyl attached to secondary carbon"

    return False, "No secondary alcohol group detected"