"""
Classifies: CHEBI:16366 anthocyanidin cation
"""
"""
Classifies: CHEBI:35218 anthocyanidin cation
"""
from rdkit import Chem

def is_anthocyanidin_cation(smiles: str):
    """
    Determines if a molecule is an anthocyanidin cation based on its SMILES string.
    Anthocyanidin cations are oxygenated derivatives of flavylium (2-phenylchromenylium).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an anthocyanidin cation, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for flavylium core (2-phenylchromenylium)
    # Updated SMARTS pattern to match fused benzene and pyrylium with phenyl substituent
    flavylium_pattern = Chem.MolFromSmarts("[o+]1c2ccccc2ccc1-c3ccccc3")
    if not mol.HasSubstructMatch(flavylium_pattern):
        return False, "Flavylium core (2-phenylchromenylium) not found"
    
    # Verify oxygen substituents (at least one oxygen atom besides the core O+)
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if oxygen_count < 2:
        return False, "Insufficient oxygen substituents (must have ≥1 beyond O+ core)"
    
    return True, "Contains flavylium core with phenyl group and oxygen substituents"