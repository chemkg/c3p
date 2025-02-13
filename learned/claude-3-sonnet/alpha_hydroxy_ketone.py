"""
Classifies: CHEBI:139588 alpha-hydroxy ketone
"""
"""
Classifies: CHEBI:18243 alpha-hydroxy ketone
An alpha-hydroxy ketone is a ketone containing a hydroxy group on the alpha-carbon relative to the C=O group.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_alpha_hydroxy_ketone(smiles: str):
    """
    Determines if a molecule is an alpha-hydroxy ketone based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alpha-hydroxy ketone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Find ketone atoms
    ketone_atoms = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and atom.GetIsAromatic() is False and atom.GetTotalDegree() == 2]
    
    # Check for hydroxy group on alpha carbon for each ketone
    for ketone_idx in ketone_atoms:
        alpha_carbon = mol.GetAtomWithIdx(ketone_idx).GetNeighbors()[0].GetIdx()
        alpha_carbon_atom = mol.GetAtomWithIdx(alpha_carbon)
        if alpha_carbon_atom.GetTotalNumHs() > 0:
            return True, "Contains a ketone with a hydroxy group on the alpha carbon"
    
    return False, "Does not contain an alpha-hydroxy ketone group"