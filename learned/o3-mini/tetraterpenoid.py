"""
Classifies: CHEBI:26935 tetraterpenoid
"""
"""
Classifies: Tetraterpenoid 
Defined as any terpenoid derived from a tetraterpene (typically with a C40 core, possibly rearranged or slightly modified),
which commonly exhibits an extended conjugated polyene chain.
Improved criteria compared to the previous version.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_tetraterpenoid(smiles: str):
    """
    Determines if a molecule is a tetraterpenoid based on its SMILES string.
    The improved classification is based on several criteria:
      1. The total number of carbon atoms must be in a narrow range (30 to 50) consistent with a C40-derived scaffold.
      2. The molecule must contain a sufficiently long conjugated chain, approximated by counting 
         non‐aromatic carbon–carbon double bonds (requiring at least 7).
      3. The molecular weight should be within a range (300 to 750 Da) typical for tetraterpenoids.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a tetraterpenoid, False otherwise.
        str: An explanation of the classification decision.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count carbon atoms
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if not (30 <= carbon_count <= 50):
        return False, f"Carbon count ({carbon_count}) is not in the expected range for tetraterpenoids (30-50)."

    # Count non-aromatic carbon-carbon double bonds
    double_bond_count = 0
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.BondType.DOUBLE:
            a1 = bond.GetBeginAtom()
            a2 = bond.GetEndAtom()
            # Only count if both atoms are carbons and are not aromatic.
            if (a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 6 and 
                not a1.GetIsAromatic() and not a2.GetIsAromatic()):
                double_bond_count += 1
    if double_bond_count < 7:
        return False, f"Found only {double_bond_count} non-aromatic C=C bonds; expected at least 7 for an extended conjugated system."

    # Calculate molecular weight
    mw = rdMolDescriptors.CalcExactMolWt(mol)
    if mw < 300 or mw > 750:
        return False, f"Molecular weight ({mw:.1f}) is not in the expected range (300-750 Da) for tetraterpenoids."

    return True, (f"Carbon count ({carbon_count}), {double_bond_count} non-aromatic C=C bonds, "
                  f"and molecular weight ({mw:.1f}) are consistent with a tetraterpenoid.")

# Example usage (uncomment the lines below to test):
# test_smiles = "CC(\\C=C\\C=C(C)C=C\\C1=C(C)C[C@@H](O)CC1(C)C)=C/C=C/C=C(C)/C=C/C=C(C)/C=C/C1=C(C)C[C@@H](O)CC1(C)C"
# print(is_tetraterpenoid(test_smiles))