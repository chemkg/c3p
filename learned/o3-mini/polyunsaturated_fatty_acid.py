"""
Classifies: CHEBI:26208 polyunsaturated fatty acid
"""
"""
Classifies: Polyunsaturated fatty acid (PUFA)
Definition: Any fatty acid containing more than one double bond.
Acids in this group are reported to have cardioprotective effects; and levels are lowered in chronic fatigue syndrome.
"""

from rdkit import Chem

def is_polyunsaturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a polyunsaturated fatty acid based on its SMILES string.
    A polyunsaturated fatty acid (PUFA) is defined as any fatty acid that contains
    more than one carbon-carbon double bond and has a carboxylic acid group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polyunsaturated fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES into an RDKit molecule object.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for the presence of a carboxylic acid group.
    # This pattern looks for a carbonyl (C=O) attached to an -OH (O).
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No carboxyl group found"
        
    # Count the number of carbon-carbon double bonds.
    # We count only C=C bonds and exclude the C=O bond from the carboxyl group.
    double_bond_count = 0
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
            atom1 = bond.GetBeginAtom()
            atom2 = bond.GetEndAtom()
            if atom1.GetAtomicNum() == 6 and atom2.GetAtomicNum() == 6:
                double_bond_count += 1
    
    # Check if we have the required minimum number of double bonds.
    if double_bond_count < 2:
        return False, f"Only {double_bond_count} carbon-carbon double bond(s) found; need at least 2"
    
    return True, f"Contains carboxyl group and {double_bond_count} carbon-carbon double bonds, meets PUFA criteria"