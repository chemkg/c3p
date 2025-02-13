"""
Classifies: CHEBI:59549 essential fatty acid
"""
"""
Classifies: Essential Fatty Acid
Definition: Any member of the sub-set of polyunsaturated fatty acid for which there is an absolute dietary requirement.
This function uses a simple approach:
    - It checks that the SMILES string can be parsed.
    - It checks for a carboxylic acid function (as free acid or as part of an ester) since fatty acids usually contain a terminal –C(=O)O group.
    - It counts the total number of carbon atoms (to ensure a "long" chain) and the number of carbon–carbon double bonds (to ensure polyunsaturation).
Note: This method is a simplified heuristic and may not capture all nuances of essential fatty acids.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_essential_fatty_acid(smiles: str):
    """
    Determines if a molecule is an essential fatty acid based on its SMILES string.
    The criteria used include:
      1. Valid molecule.
      2. Presence of a carboxylic acid or ester carboxyl group.
      3. A minimum number of carbon atoms in the molecule (indicative of a long hydrocarbon chain).
      4. At least two carbon–carbon double bonds indicating polyunsaturation.
      5. A molecular weight consistent with a fatty acid (not too low).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as an essential fatty acid, False otherwise.
        str: Reason for the classification decision.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # Check for the carboxylic acid functionality.
    # The SMARTS pattern here captures a carbonyl bonded to an oxygen that is either acidic (HO) or has a negative charge.
    acid_pattern = Chem.MolFromSmarts("C(=O)[O;H,-]")
    if not mol.HasSubstructMatch(acid_pattern):
        return False, "No carboxylic acid (or related ester acid) functional group found."
    
    # Count total carbon atoms in the molecule.
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 8:
        return False, f"Too few carbons ({c_count}) to be a fatty acid chain."
    
    # Count the number of carbon-carbon double bonds.
    double_bond_count = 0
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.BondType.DOUBLE:
            a1 = bond.GetBeginAtom()
            a2 = bond.GetEndAtom()
            if a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 6:
                double_bond_count += 1
    if double_bond_count < 2:
        return False, f"Insufficient polyunsaturation (found {double_bond_count} C=C bonds, need at least 2)."
    
    # Check molecular weight (a very short chain would be too light to be a fatty acid).
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 150:
        return False, f"Molecular weight is too low ({mol_wt:.1f} Da) for a fatty acid."
    
    return True, "Molecule contains a carboxylic acid group, has a long hydrocarbon chain, and is polyunsaturated."