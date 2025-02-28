"""
Classifies: CHEBI:65111 3-substituted propionyl-CoA(4-)
"""
"""
Classifies: 3-substituted propionyl-CoA(4-)
"""
from rdkit import Chem

def is_3_substituted_propionyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is a 3-substituted propionyl-CoA(4-) based on its SMILES string.
    This involves finding the Coenzyme A structure and ensuring substitution on the third carbon
    of the propionyl group.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-substituted propionyl-CoA(4-), False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the general CoA backbone SMARTS pattern, updated for substitution on the third carbon
    coa_pattern = Chem.MolFromSmarts("C(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No CoA backbone found"

    # Define the 3-substituted propionyl pattern specifically on the third carbon
    # Check for a generic R group connected to the third carbon atom
    propionyl_substituted_pattern = Chem.MolFromSmarts("C(=O)SCC([#6,#7,#8,#9])")
    if not mol.HasSubstructMatch(propionyl_substituted_pattern):
        return False, "No 3-substituted propionyl group found on the third position"

    # If both patterns match, assume it is a 3-substituted propionyl-CoA(4-)
    return True, "Contains CoA backbone with a 3-substituted propionyl group on the third carbon"