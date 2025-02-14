"""
Classifies: CHEBI:29256 thiol
"""
"""
Classifies: CHEBI:29232 thiol
"""
from rdkit import Chem

def is_thiol(smiles: str):
    """
    Determines if a molecule is a thiol based on its SMILES string.
    A thiol is an organosulfur compound in which a thiol group (-SH) is attached to a carbon atom
    of any aliphatic or aromatic moiety.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a thiol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens
    mol = Chem.AddHs(mol)
    
    # Define thiol SMARTS pattern: sulfur bonded to hydrogen and carbon
    thiol_pattern = Chem.MolFromSmarts("[#6]-[SH]")
    if thiol_pattern is None:
        return False, "Invalid SMARTS pattern for thiol"
    
    # Search for thiol group
    if mol.HasSubstructMatch(thiol_pattern):
        return True, "Contains a thiol group (-SH) attached to a carbon atom"
    else:
        return False, "No thiol group (-SH) attached to a carbon atom found"