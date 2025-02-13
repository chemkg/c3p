"""
Classifies: CHEBI:25029 leukotriene
"""
#!/usr/bin/env python
"""
Classifies: leukotriene (icosanoids derived from arachidonic acid with 4 double bonds,
   of which 3 are conjugated)
"""

from rdkit import Chem

def is_leukotriene(smiles: str):
    """
    Determines if a molecule is a leukotriene based on its SMILES string.
    A leukotriene is defined here as a C20-derived (icosanoid) molecule originating from
    arachidonic acid and having exactly 4 carbon–carbon double bonds, with a contiguous
    system of three conjugated double bonds.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule meets the leukotriene criteria, False otherwise.
        str: Explanation of the classification.
    """
    # Parse SMILES string into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check that the molecule has at least 20 carbon atoms as expected for C20 icosanoids
    carbon_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    if len(carbon_atoms) < 20:
        return False, f"Only {len(carbon_atoms)} carbon atoms found; expected at least 20 for an icosanoid"

    # Count carbon-carbon double bonds.
    # We only count bonds where both atoms are carbons.
    c_c_double_bond_count = 0
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
            a1, a2 = bond.GetBeginAtom(), bond.GetEndAtom()
            if a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 6:
                c_c_double_bond_count += 1
    if c_c_double_bond_count != 4:
        return False, f"Found {c_c_double_bond_count} carbon–carbon double bonds; expected exactly 4"
    
    # Search for a contiguous system of three conjugated double bonds.
    # This SMARTS pattern looks for three double bonds separated by single bonds: C=C-C=C-C=C.
    conjugated_pattern = Chem.MolFromSmarts("[C]=[C]-[C]=[C]-[C]=[C]")
    if not mol.HasSubstructMatch(conjugated_pattern):
        return False, "No contiguous tri-conjugated diene system (three conjugated C=C bonds) found"
    
    return True, "Molecule has a C20 or greater backbone with exactly 4 C=C bonds, including a contiguous system of 3 conjugated double bonds typical of leukotrienes"


# For quick testing, you could uncomment the lines below:
# test_smiles = "CCCC/C=C/C[C@@H](O)/C=C/C=C/C=C/[C@@H](O)CCCC(O)=O"  # 6-trans-leukotriene B4
# print(is_leukotriene(test_smiles))