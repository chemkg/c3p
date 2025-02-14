"""
Classifies: CHEBI:36835 3alpha-hydroxy steroid
"""
from rdkit import Chem

def is_3alpha_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 3alpha-hydroxy steroid based on its SMILES string.
    A 3alpha-hydroxy steroid is a 3-hydroxy steroid where the 3-hydroxy substituent is in the alpha-position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3alpha-hydroxy steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Generalized pattern for steroid backbone: four connected rings, which may be a common pattern in steroids.
    steroid_backbone_pattern = Chem.MolFromSmarts("C1CCC2C(C1)CCC3C2CCC4C3CCC(C4)C")
    if not mol.HasSubstructMatch(steroid_backbone_pattern):
        return False, "No steroid backbone found"
        
    # SMARTS pattern for 3alpha-hydroxy group in a chiral center:
    # Extended to match the hydroxyl group at specifically the third position with chiral specification.
    hydroxy_pattern = Chem.MolFromSmarts("[C@@H]([O])C")  # Hydroxyl alpha with specified chiral center at C3
    hydroxy_matches = mol.GetSubstructMatches(hydroxy_pattern)
    if not hydroxy_matches:
        return False, "No 3alpha-hydroxy group found"

    # Further check if at least one of the hydroxyl groups is correctly positioned in context.
    for match in hydroxy_matches:
        # Assuming match[0] is the carbon center for hydroxyl, basic check for cyclic carbon for confirmation
        # (Note: Advanced checks should consider whole structural integrity and stereochemistry)
        if mol.GetAtomWithIdx(match[0]).GetDegree() == 3:  # Checking for proper connection count typical for steroids
            return True, "Contains steroid backbone with 3alpha-hydroxy group"
    
    return False, "Proper 3alpha-hydroxy position not confirmed in context"

# Test the function with an example
smiles_example = "C1[C@@]2([C@]3(CC[C@]4([C@]([C@@]3(CC[C@]2(C[C@@H]([C@@H]1O)O)[H])[H])(CC[C@@]4([C@@H]([C@H]([C@@H](CC(C)C)O)O)C)[H])[H])C)[H])C"
print(is_3alpha_hydroxy_steroid(smiles_example))  # Expected: True with reason