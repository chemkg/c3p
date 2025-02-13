"""
Classifies: CHEBI:26935 tetraterpenoid
"""
"""
Classifies: Tetraterpenoid 
Defined as any terpenoid derived from a tetraterpene (typically with a C40 core, possibly rearranged or slightly modified),
which commonly exhibits an extended conjugated polyene chain.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_tetraterpenoid(smiles: str):
    """
    Determines if a molecule is a tetraterpenoid based on its SMILES string.
    The classification is based on several criteria:
      1. The total number of carbon atoms should be in a range consistent with a C40
         skeleton – we allow some flexibility (here 30 to 50 carbons).
      2. The molecule should possess several non‐aromatic carbon–carbon double bonds,
         as a proxy for an extended conjugated system (typical in carotenoid derivatives).
      3. A minimum molecular weight is expected.
    
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

    # Count carbon atoms in the molecule.
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    
    # Check if carbon count is roughly consistent with a tetraterpene origin.
    # Tetraterpenoids are derived from eight isoprene units (8 x 5 = 40 carbons)
    # but modifications can cause slight deviations.
    if not (30 <= carbon_count <= 50):
        return False, f"Carbon count ({carbon_count}) is not in the expected range for tetraterpenoids (30-50)."
    
    # Count non-aromatic carbon-carbon double bonds as an (approximate) measure for an extended conjugated system.
    double_bond_count = 0
    for bond in mol.GetBonds():
        # Check if bond is a double bond and connects two carbons
        if bond.GetBondType() == Chem.BondType.DOUBLE:
            a1 = bond.GetBeginAtom()
            a2 = bond.GetEndAtom()
            if a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 6:
                # We ignore bonds if they involve aromatic systems (since conjugated polyenes in tetraterpenoids are non-aromatic)
                if not a1.GetIsAromatic() and not a2.GetIsAromatic():
                    double_bond_count += 1
    
    # Many tetraterpenoids (such as carotenoids) display an extensive conjugated system.
    # While the exact number can vary, we expect several double bonds.
    if double_bond_count < 5:
        return False, f"Found only {double_bond_count} non-aromatic C=C bonds; expected an extended conjugated system."
    
    # Check the molecular weight as a further rough indicator.
    mw = rdMolDescriptors.CalcExactMolWt(mol)
    if mw < 300:
        return False, f"Molecular weight {mw:.1f} is too low for a typical tetraterpenoid."
    
    # Passed basic criteria for a tetraterpenoid candidate
    return True, f"Carbon count ({carbon_count}), {double_bond_count} C=C bonds, and molecular weight ({mw:.1f}) are consistent with a tetraterpenoid."
    
# Example usage (uncomment below lines to test):
# test_smiles = "CC(\C=C\C=C(C)\C=C\C1=C(C)C[C@@H](O)CC1(C)C)=C/C=C/C=C(C)/C=C/C=C(C)/C=C/C1=C(C)C[C@@H](O)CC1(C)C"  # zeaxanthin example
# print(is_tetraterpenoid(test_smiles))