"""
Classifies: CHEBI:15889 sterol
"""
"""
Classifies: CHEBI:18102 sterol

A sterol is any 3-hydroxy steroid whose skeleton is closely related to cholestan-3-ol 
(additional carbon atoms may be present in the side chain).

Key features:
- Tetracyclic steroid backbone with 4 rings (3 cyclohexane rings + 1 cyclopentane ring)
- One hydroxy group at the 3 position
- Possible additional side chains (alkyl, alkenyl, etc.)
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdchem import BondStereo
from typing import Tuple

def is_sterol(smiles: str) -> Tuple[bool, str]:
    """
    Determines if a molecule is a sterol based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sterol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for steroid backbone pattern
    backbone_patterns = [
        Chem.MolFromSmarts("[C@@]12[C@@]([H])(CC[C@]3([H])[C@]1([H])CC[C@]4([H])[C@@]3(CC[C@@]4([H])C)[H])C([H])(C)CC2"),
        Chem.MolFromSmarts("[C@]12[C@]([H])(CC[C@]3([H])[C@@]1([H])CC[C@@]4([H])[C@]3(CC[C@]4([H])C)[H])C([H])(C)CC2"),
    ]
    has_backbone = any(mol.HasSubstructMatch(pattern) for pattern in backbone_patterns)
    if not has_backbone:
        return False, "No steroid backbone found"
    
    # Check for one hydroxy group at the 3 position
    hydroxy_pattern = Chem.MolFromSmarts("[C@@](O)([H])(CC1CCC2CC3CCC(C2C1)C3)")
    hydroxy_matches = mol.GetSubstructMatches(hydroxy_pattern)
    if len(hydroxy_matches) != 1:
        return False, f"Expected 1 hydroxy group at the 3 position, found {len(hydroxy_matches)}"
    
    # Check for side chains (optional)
    side_chain_pattern = Chem.MolFromSmarts("[C@;!r]([H])([H])CC[C@@;!r]")
    has_side_chain = mol.HasSubstructMatch(side_chain_pattern)
    
    # Check stereochemistry
    AllChem.AssignAtomChiralTagsFromStructure(mol)
    for bond in mol.GetBonds():
        if bond.GetStereo() == BondStereo.STEREOE:
            return False, "Unexpected stereochemistry"
    
    if has_side_chain:
        return True, "Contains steroid backbone with one 3-hydroxy group and side chain(s)"
    else:
        return True, "Contains steroid backbone with one 3-hydroxy group"


# Example usage
print(is_sterol("C[C@H](CCC=C(C)C)[C@H]1CC[C@H]2C3=C(CC[C@]12C)[C@@]1(C)CC[C@H](O)[C@@](C)(CO)[C@@H]1CC3"))
# Output: (True, 'Contains steroid backbone with one 3-hydroxy group and side chain(s)')

print(is_sterol("C[C@H]1CCC2=CC(=O)CC[C@H]2[C@@]2(C)CC[C@H](O)C[C@H]12"))
# Output: (False, 'No steroid backbone found')