"""
Classifies: CHEBI:18310 alkane
"""
#!/usr/bin/env python
"""
Classifies: alkane (acyclic saturated hydrocarbon having formula CnH2n+2)
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_alkane(smiles: str):
    """
    Determines if a molecule is an alkane based on its SMILES string.
    An alkane is defined as an acyclic (non-cyclic) hydrocarbon that contains
    only carbon (C) and hydrogen (H) atoms, has only single bonds and satisfies
    the formula CnH2n+2.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an alkane, False otherwise
        str: Reason for classification
    """
    # Parse SMILES using RDKit
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check that the molecule consist only of carbon and hydrogen atoms
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num not in (1, 6):  # Only H (1) and C (6) are allowed
            return False, f"Contains element with atomic number {atomic_num}, not allowed in an alkane"

    # Check that the molecule is acyclic (no rings)
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() > 0:
        return False, "Molecule contains rings, hence not an acyclic alkane"

    # Check that all bonds are single bonds and non-aromatic 
    for bond in mol.GetBonds():
        # BondType.SINGLE represents a single bond; also check if the bond is marked aromatic
        if bond.GetBondType() != Chem.rdchem.BondType.SINGLE or bond.GetIsAromatic():
            return False, "Molecule contains bonds that are not single (or are aromatic)"

    # Add explicit hydrogens to accurately count them
    mol_with_H = Chem.AddHs(mol)

    # Count carbons and hydrogens
    c_count = sum(1 for atom in mol_with_H.GetAtoms() if atom.GetAtomicNum() == 6)
    h_count = sum(1 for atom in mol_with_H.GetAtoms() if atom.GetAtomicNum() == 1)

    # Check if the formula matches CnH2n+2
    if h_count != 2 * c_count + 2:
        return False, f"Molecular formula does not match CnH2n+2 (found C{c_count}H{h_count})"
    
    return True, "Molecule is an acyclic saturated hydrocarbon that matches CnH2n+2 (alkane)"

# Example tests (uncomment to run)
# test_smiles = "CCCCCC"  # hexane
# result, reason = is_alkane(test_smiles)
# print(result, reason)