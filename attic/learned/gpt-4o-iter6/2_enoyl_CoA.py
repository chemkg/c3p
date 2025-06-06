"""
Classifies: CHEBI:19573 2-enoyl-CoA
"""
from rdkit import Chem


def is_2_enoyl_CoA(smiles: str):
    """
    Determines if a molecule is a 2-enoyl-CoA based on its SMILES string.
    A 2-enoyl-CoA is defined as an unsaturated fatty acyl-CoA in which the S-acyl group contains a double bond between positions 2 and 3.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2-enoyl-CoA, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the pattern for Coenzyme A with a 2-enoyl feature
    # Thioester linkage with a double bond on the next fatty acyl chain (C=C-C(=O) or similar)
    
    # Extremely simplified pattern for 2-enoyl thioester linkage
    enoyl_pattern = Chem.MolFromSmarts("C(=O)SC=CC")
    
    # Check if molecule contains the 2-enoyl feature structure
    if not mol.HasSubstructMatch(enoyl_pattern):
        return False, "No 2-enoyl pattern found"
        
    # Ensure CoA structure is present
    # This is checking for a more general context of CoA structure
    coA_pattern = Chem.MolFromSmarts("COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H](O[C@H]([C@@H]1OP(O)(O)=O)n2cnc3c(N)ncnc23)")
    if not mol.HasSubstructMatch(coA_pattern):
        return False, "No Coenzyme A backbone detected"
    
    return True, "Contains 2-enoyl feature with Coenzyme A structure"