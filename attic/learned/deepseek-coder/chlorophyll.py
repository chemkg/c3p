"""
Classifies: CHEBI:28966 chlorophyll
"""
"""
Classifies: CHEBI:35255 chlorophyll
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_chlorophyll(smiles: str):
    """
    Determines if a molecule is a chlorophyll based on its SMILES string.
    A chlorophyll is a magnesium porphyrin with a fifth ring and typically a long phytol chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a chlorophyll, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for magnesium in the molecule
    mg_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 12]
    if len(mg_atoms) != 1:
        return False, "No magnesium atom found"

    # Define the porphyrin core pattern with the fifth ring
    porphyrin_pattern = Chem.MolFromSmarts("[Mg]1234n1c(c2)c(c3)c(c4)c5c1c(c2)c(c3)c(c4)c5")
    if not mol.HasSubstructMatch(porphyrin_pattern):
        return False, "No porphyrin core with fifth ring found"

    # Check for a long phytol chain (at least 10 carbons)
    phytol_pattern = Chem.MolFromSmarts("[CX4][CX4][CX4][CX4][CX4][CX4][CX4][CX4][CX4][CX4]")
    if not mol.HasSubstructMatch(phytol_pattern):
        return False, "No long phytol chain found"

    # Check molecular weight - chlorophylls typically >800 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 800:
        return False, "Molecular weight too low for chlorophyll"

    return True, "Contains magnesium porphyrin core with fifth ring and long phytol chain"