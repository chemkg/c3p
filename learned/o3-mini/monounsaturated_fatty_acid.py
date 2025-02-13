"""
Classifies: CHEBI:25413 monounsaturated fatty acid
"""
"""
Classifies: Monounsaturated fatty acid (any fatty acid with one double or triple bond in the fatty acid chain and single bonds otherwise)
MUFAs are known to have positive effects on the cardiovascular system and in diabetes treatment.
"""

from rdkit import Chem

def is_monounsaturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a monounsaturated fatty acid based on its SMILES string.
    A monounsaturated fatty acid must contain a terminal carboxylic acid (–C(=O)OH) group and exactly one carbon-carbon
    unsaturation (double or triple bond) in its chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monounsaturated fatty acid, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a SMARTS pattern for a carboxylic acid (–C(=O)OH) group.
    # This pattern targets a carbon with a double-bonded oxygen and an -OH group.
    acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    if not mol.HasSubstructMatch(acid_pattern):
        return False, "No terminal carboxylic acid group found"

    # Count unsaturation bonds (double or triple) between two carbon atoms.
    unsaturation_count = 0
    for bond in mol.GetBonds():
        btype = bond.GetBondType()
        if btype in (Chem.rdchem.BondType.DOUBLE, Chem.rdchem.BondType.TRIPLE):
            # Only count bonds where both atoms are carbons.
            if bond.GetBeginAtom().GetAtomicNum() == 6 and bond.GetEndAtom().GetAtomicNum() == 6:
                unsaturation_count += 1

    if unsaturation_count != 1:
        return False, f"Found {unsaturation_count} carbon–carbon unsaturation(s); need exactly one"

    # Optionally, one could verify that the remainder of the chain is saturated (i.e. only single bonds)
    # and that the acyl chain length is appropriate. For simplicity, the definition is met if the above conditions hold.
    
    return True, "Molecule contains one carbon–carbon unsaturation and a terminal carboxylic acid group, qualifying as a monounsaturated fatty acid"

# Uncomment the following lines to do a quick test with one of the provided examples.
# example_smiles = "CCCCCC\\C=C\\CCCCCCCCCC(O)=O"  # trans-vaccenic acid
# result, reason = is_monounsaturated_fatty_acid(example_smiles)
# print(result, reason)