"""
Classifies: CHEBI:73155 trienoic fatty acid
"""
"""
Classifies: Trienoic fatty acid
Definition: Any polyunsaturated fatty acid (only C, H, O) that contains three carbon–carbon double bonds 
and has a free terminal carboxylic acid group.
Heuristics used:
  1. Only atoms allowed are C, H, and O.
  2. The molecule must contain a free (i.e. non‐ester) carboxylic acid group (matched by SMARTS for both protonated
     and deprotonated forms). In that match the acid carbon is required to be terminal (only one carbon neighbor).
  3. The molecule must contain exactly three C=C double bonds (ignoring any carbonyl bonds).
"""

from rdkit import Chem

def is_trienoic_fatty_acid(smiles: str):
    """
    Determines if a molecule is a trienoic fatty acid based on its SMILES string.
    A trienoic fatty acid (by our heuristics) must:
      • contain only C, H, and O atoms,
      • possess a free, terminal carboxylic acid group, and
      • have exactly three carbon–carbon double bonds.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule meets our criteria, False otherwise.
        str: A reason message explaining the outcome.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Reject molecules containing atoms other than C, H, and O.
    allowed_atomic_nums = {1, 6, 8}
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in allowed_atomic_nums:
            return False, "Contains atoms other than C, H, and O; not a fatty acid"

    # Look for a free (non-ester) carboxylic acid group.
    # We use two SMARTS variants to include protonated and deprotonated forms.
    acid_smarts_list = [
        "[CX3](=O)[OX2H1]",  # protonated acid, e.g. C(=O)O
        "[CX3](=O)[O-]"      # deprotonated acid, e.g. C(=O)[O-]
    ]
    acid_match = None
    for smarts in acid_smarts_list:
        acid_pat = Chem.MolFromSmarts(smarts)
        matches = mol.GetSubstructMatches(acid_pat)
        if matches:
            # Use the first match.
            acid_match = matches[0]
            break
    if acid_match is None:
        return False, "No carboxylic acid group found; not a fatty acid"
    # By convention in our SMARTS the first atom is the acid (carbonyl) carbon.
    acidC_idx = acid_match[0]
    acidC = mol.GetAtomWithIdx(acidC_idx)
    # Check that the acid carbon is terminal (i.e. attached to exactly one carbon atom).
    carbon_neighbors = [nbr for nbr in acidC.GetNeighbors() if nbr.GetAtomicNum() == 6]
    if len(carbon_neighbors) != 1:
        return False, "Carboxylic acid group is not terminal; may be esterified or not a free fatty acid"
    # Also check that the alpha carbon (neighboring carbon) is not overly substituted.
    alphaC = carbon_neighbors[0]
    # Count explicit hydrogens (use GetTotalNumHs which returns total of implicit and explicit)
    alphaH = alphaC.GetTotalNumHs()
    # In a simple fatty acid, the alpha carbon is usually CH2 so we expect 2 hydrogens.
    if alphaH != 2:
        return False, "Alpha carbon adjacent to acid group is over/substituted; not a typical fatty acid backbone"
    
    # Count the number of pure carbon–carbon double bonds.
    double_bond_count = 0
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
            a1 = bond.GetBeginAtom()
            a2 = bond.GetEndAtom()
            # Only count bond if both atoms are carbon.
            if a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 6:
                double_bond_count += 1

    if double_bond_count != 3:
        return False, f"Found {double_bond_count} carbon–carbon double bonds; need exactly 3"
    
    # Passed all tests.
    return True, "Molecule is a trienoic fatty acid (free terminal carboxylic acid group and exactly three C=C bonds)"

# Example testing (you can remove or comment these out in production)
if __name__ == "__main__":
    test_smiles = [
        "OC(=O)CCCCCCC/C=C\\CC(=O)/C=C(\\O)/C=C\\CCO",  # True positive example
        "CCCC\\C=C/C\\C=C/CC(O)C(O)C\\C=C/CCCC([O-])=O"  # False positive example expected to fail our extra tests
    ]
    for smi in test_smiles:
        result, reason = is_trienoic_fatty_acid(smi)
        print(f"SMILES: {smi}\n  Result: {result}\n  Reason: {reason}\n")