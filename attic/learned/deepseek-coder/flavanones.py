"""
Classifies: CHEBI:28863 flavanones
"""
"""
Classifies: CHEBI:72564 flavanone
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_flavanone(smiles: str):
    """
    Determines if a molecule is a flavanone based on its SMILES string.
    A flavanone has a 3,4-dihydro-2-aryl-2H-1-benzopyran-4-one skeleton.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a flavanone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the core flavanone pattern
    flavanone_pattern = Chem.MolFromSmarts("[C@@H]1CC(=O)c2ccccc2O1")
    if not mol.HasSubstructMatch(flavanone_pattern):
        return False, "No flavanone core structure found"

    # Check for aromatic ring at position 2
    aryl_pattern = Chem.MolFromSmarts("[C@@H]1CC(=O)c2ccccc2O1-c3ccccc3")
    if not mol.HasSubstructMatch(aryl_pattern):
        return False, "No aryl group at position 2"

    # Verify the 3,4-dihydro structure (no double bond between C3 and C4)
    dihydro_pattern = Chem.MolFromSmarts("[C@H]1CC(=O)c2ccccc2O1")
    if not mol.HasSubstructMatch(dihydro_pattern):
        return False, "No 3,4-dihydro structure found"

    # Count oxygen atoms (should be at least 2: one in the ring, one in the ketone)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 2:
        return False, "Insufficient oxygen atoms"

    # Count carbon atoms (should be at least 15 for basic flavanone structure)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 15:
        return False, "Insufficient carbon atoms"

    # Check molecular weight (should be >200 for flavanones)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200:
        return False, "Molecular weight too low for flavanone"

    return True, "Contains 3,4-dihydro-2-aryl-2H-1-benzopyran-4-one skeleton"