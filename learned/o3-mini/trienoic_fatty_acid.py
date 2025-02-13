"""
Classifies: CHEBI:73155 trienoic fatty acid
"""
"""
Classifies: Trienoic fatty acid 
Definition: Any polyunsaturated fatty acid that contains three double bonds.
A candidate must have a carboxylic acid group and exactly three C=C double bonds (non-carbonyl).
"""
from rdkit import Chem

def is_trienoic_fatty_acid(smiles: str):
    """
    Determines if a molecule is a trienoic fatty acid based on its SMILES string.
    A trienoic fatty acid must be a fatty acid (contain a carboxylic acid group)
    and have three carbon–carbon double bonds.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if molecule is a trienoic fatty acid, False otherwise.
        str: Reason for the classification.
    """
    # Parse SMILES into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS for a carboxylic acid group.
    # This pattern matches both protonated carboxylic acids and deprotonated forms.
    acid_smarts_list = [
        "[CX3](=O)[OX2H1]",  # protonated acid, e.g., C(=O)O
        "[CX3](=O)[O-]"      # deprotonated acid, e.g., C(=O)[O-]
    ]
    has_acid = False
    for smarts in acid_smarts_list:
        acid_pattern = Chem.MolFromSmarts(smarts)
        if mol.HasSubstructMatch(acid_pattern):
            has_acid = True
            break
    if not has_acid:
        return False, "No carboxylic acid group found; not a fatty acid"

    # Count the number of C=C double bonds excluding those in a carbonyl group.
    # We consider only double bonds between two carbons.
    double_bond_count = 0
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
            a1 = bond.GetBeginAtom()
            a2 = bond.GetEndAtom()
            # Check if both atoms are carbons
            if a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 6:
                double_bond_count += 1

    if double_bond_count != 3:
        return False, f"Found {double_bond_count} carbon–carbon double bonds; need exactly 3"

    # Optionally, one may check overall chain length or other descriptors,
    # but based on the definition given, the above checks suffice.
    return True, "Molecule is a trienoic fatty acid (contains a carboxylic acid group and three C=C bonds)"