"""
Classifies: CHEBI:74716 withanolide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_withanolide(smiles: str):
    """
    Determines if a molecule is a withanolide based on its SMILES string.
    A withanolide is a C28 steroid with a modified side chain forming a lactone ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a withanolide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Detect the steroid core using a more specific pattern (tetracyclic ring system).
    steroid_pattern = Chem.MolFromSmarts("[C]1[C][C]2[C]([C]1)[C]3[C]([C]2)[C][C]4[C]([C]3)[C][C]([C]4)[C]")
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid core found"
    
    # 2. Detect the lactone ring
    lactone_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    if not mol.HasSubstructMatch(lactone_pattern):
         return False, "No lactone ring found"
    
    # 3. Verify C28 steroid (approximately 28 carbons)
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 26 or carbon_count > 30:
        return False, f"Carbon count ({carbon_count}) not in range for a withanolide (26-30)"
    
    # 4. Check for a carbonyl in sidechain.
    carbonyl_pattern = Chem.MolFromSmarts("[CX3]=[OX1]")
    
    # 5. Check for the side chain attached to core with a carbonyl.
    side_chain_pattern = Chem.MolFromSmarts("[C]~[C]~[C]~[CX3]=[OX1]")
    if not mol.HasSubstructMatch(side_chain_pattern):
        side_chain_pattern = Chem.MolFromSmarts("[C]~[C]~[C]~[C]~[CX3]=[OX1]")
        if not mol.HasSubstructMatch(side_chain_pattern):
            return False, "No characteristic side chain with carbonyl found"
    
    return True, "Contains a steroid core with lactone and a characteristic sidechain"