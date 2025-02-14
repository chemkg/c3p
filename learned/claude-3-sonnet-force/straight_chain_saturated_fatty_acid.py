"""
Classifies: CHEBI:39418 straight-chain saturated fatty acid
"""
"""
Classifies: CHEBI:35838 straight-chain saturated fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, SanitizeMol

def is_straight_chain_saturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a straight-chain saturated fatty acid based on its SMILES string.
    A straight-chain saturated fatty acid is a saturated fatty acid lacking a side-chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a straight-chain saturated fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Remove explicit hydrogens and parentheses for normalization
    SanitizeMol(mol, sanitizeOps=Chem.SanitizeFlags.SANITIZE_PROPERTIES |
                             Chem.SanitizeFlags.SANITIZE_FINDRADICALS |
                             Chem.SanitizeFlags.SANITIZE_KEKULIZE |
                             Chem.SanitizeFlags.SANITIZE_SETAROMATICITY)
    Chem.RemoveHs(mol)
    Chem.RemoveParentheses(mol)
    
    # Check for carboxylic acid group (-C(=O)O)
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)[O;H,-]")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No carboxylic acid group found"
    
    # Check for straight carbon chain (no branches)
    straight_chain_pattern = Chem.MolFromSmarts("C(C)(C)C")
    if mol.HasSubstructMatch(straight_chain_pattern):
        return False, "Branched carbon chain found"
    
    # Check for saturation (no double or triple bonds)
    for bond in mol.GetBonds():
        if bond.GetBondType() not in (Chem.BondType.SINGLE, Chem.BondType.AROMATIC):
            if bond.GetIsotopeSum() == 0:  # Ignore isotopic labeling
                return False, "Unsaturated bond found"
    
    # Check for hydroxy groups (optional, but no other substituents)
    hydroxy_pattern = Chem.MolFromSmarts("O[H]")
    other_substituents = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in (1, 6, 8):  # H, C, O
            other_substituents.append(atom.GetSymbol())
    if other_substituents:
        return False, f"Found other substituents: {', '.join(other_substituents)}"
    
    # Additional check for typical molecular weight or chain length
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 100 or mol_wt > 600:
        return False, "Molecular weight outside typical range for fatty acids"
    
    return True, "Straight-chain saturated fatty acid"