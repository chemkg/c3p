"""
Classifies: CHEBI:17389 2-monoglyceride
"""
"""
Classifies: 2-monoglyceride, defined as 'A monoglyceride in which the acyl substituent is located at position 2.'
The key feature is a glycerol backbone with the central carbon bearing exactly one acyl ester 
group (–O–C(=O)R) and two terminal CH2OH groups.
Improvements:
  • We refined the SMARTS to require the central carbon has exactly one hydrogen (i.e. [C;H1]) and
    that it is bound to two CH2OH groups and one ester oxygen.
  • We additionally check that the overall molecule does not contain phosphorus (to avoid phospholipids).
"""

from rdkit import Chem

def is_2_monoglyceride(smiles: str):
    """
    Determines if a molecule is a 2-monoglyceride based on its SMILES string.
    The characteristic fragment is:
         CH2OH – CH(OC(=O)R) – CH2OH
    where the central carbon is secondary. In addition, the molecule should not contain phosphorus.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is identified as a 2-monoglyceride.
        str: Explanation for the classification.
    """
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Avoid false positives: if any phosphorus atoms are present (common in phospholipids) then reject.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 15:
            return False, "Phosphorus detected; molecule is likely not a simple 2-monoglyceride"
    
    # Define a refined SMARTS pattern for the 2-monoglyceride substructure.
    # The pattern looks for a central secondary carbon ([C;H1]) that is bound to:
    #  - Two CH2OH groups: ([CH2]O)
    #  - One ester oxygen leading to an acyl group: O[C](=O)[*]
    # Note: [*] ensures that we match any atom at the end of the acyl chain.
    pattern_smarts = "[C;H1]([CH2]O)([CH2]O)O[C](=O)[*]"
    pattern = Chem.MolFromSmarts(pattern_smarts)
    if pattern is None:
        return False, "Could not compile SMARTS pattern"
    
    # Get all matches of the pattern in the molecule
    matches = mol.GetSubstructMatches(pattern)
    if not matches:
        return False, "2-monoglyceride substructure not found"
    
    # For a proper 2-monoglyceride, we expect exactly one occurrence of this fragment.
    if len(matches) != 1:
        return False, ("Multiple 2-monoglyceride-like fragments found; "
                       "the molecule may be di-/tri-acylated or not a simple 2-monoglyceride")
    
    return True, "Molecule contains a glycerol backbone with a single acyl ester at the 2-position"

# For testing purposes (can be removed in production):
if __name__ == "__main__":
    test_smiles = [
        "O=C(CCC/C=C\\C\\C=C/C=C/C(C/C=C\\CCCCC)OO)OC(CO)CO",  # 12-HPETE 2-glyceryl ester (expected: True)
        "CCCC\\C=C/C\\C=C/CCCCCCCC(=O)OC(CO)CO",                # 2-linoleoylglycerol (expected: True)
        "CCCCCCCC\\C=C/CCCCCCCC(=O)OC(CO)CO",                  # 2-oleoylglycerol (expected: True)
        "CCCCCCCCCCCCCCCC(=O)OC(CO)CO",                        # 2-palmitoylglycerol (expected: True)
        "[H]C(CO)(CO)OC(=O)CCC\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCCC",  # 2-arachidonoylglycerol (expected: True)
        # A known false positive: a glycerol backbone in a phospholipid should be rejected due to phosphorus.
        "CCCCCCCC\\C=C/CCCCCCCC(=O)OC[C@H](COP([O-])(=O)OC[C@H]([NH3+])C([O-])=O)OC(=O)CCC\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCCC"
    ]
    
    for s in test_smiles:
        flag, reason = is_2_monoglyceride(s)
        print(f"SMILES: {s}\nResult: {flag}\nReason: {reason}\n")