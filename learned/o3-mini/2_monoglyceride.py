"""
Classifies: CHEBI:17389 2-monoglyceride
"""
"""
Classifies: 2-monoglyceride, defined as 'A monoglyceride in which the acyl substituent is located at position 2.'
The key feature is a glycerol backbone with a single acyl ester at the middle carbon:
     CH2OH – CH(OC(=O)R) – CH2OH
The algorithm:
  1. Parse the SMILES string.
  2. Reject molecules containing phosphorus (to avoid phospholipids).
  3. Look for the refined SMARTS pattern for a 2-monoglyceride fragment.
  4. Confirm that the two terminal oxygen atoms in the CH2O groups are free alcohols,
     by checking that each has at least one hydrogen attached.
  5. If exactly one matching fragment is found and the additional checks pass,
     the molecule is classified as a 2-monoglyceride.
"""

from rdkit import Chem

def is_2_monoglyceride(smiles: str):
    """
    Determines if a molecule is a 2-monoglyceride based on its SMILES string.
    Characteristic fragment: CH2OH – CH(OC(=O)R) – CH2OH, with the acyl substituent on the middle carbon.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is identified as a 2-monoglyceride, False otherwise.
        str: Explanation for the classification.
    """
    # Convert the SMILES string to an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Reject molecules containing phosphorus (e.g., phospholipids).
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 15:
            return False, "Phosphorus detected; molecule is likely not a simple 2-monoglyceride"

    # Define a SMARTS pattern for the 2-monoglyceride fragment.
    # It searches for a secondary carbon ([C;H1]) attached to:
    #   - Two CH2O groups ([CH2]O) (i.e. free -CH2OH), and
    #   - One acyl ester group O[C](=O)[*].
    pattern_smarts = "[C;H1]([CH2]O)([CH2]O)O[C](=O)[*]"
    pattern = Chem.MolFromSmarts(pattern_smarts)
    if pattern is None:
        return False, "Could not compile SMARTS pattern"
    
    # Search for matches of the fragment.
    matches = mol.GetSubstructMatches(pattern)
    if not matches:
        return False, "2-monoglyceride substructure not found"
    
    # We expect a single occurrence of the fragment in a simple 2-monoglyceride.
    if len(matches) != 1:
        return False, ("Multiple 2-monoglyceride-like fragments found; " +
                       "the molecule might be di-/tri-acylated or not a simple 2-monoglyceride")
    
    # The order in our SMARTS is:
    # index 0: the central glycerol carbon
    # index 1: the first CH2O group's carbon (terminal group)
    # index 2: the second CH2O group's carbon (terminal group)
    # index 3: the oxygen from the ester bond (leading to acyl group)
    # We now verify that the two terminal CH2O groups are indeed free alcohols.
    match = matches[0]
    # For each of the CH2O groups, get the oxygen neighbor directly attached.
    # The SMARTS [CH2]O ensures that in each terminal branch a CH2 is connected to an oxygen.
    # We need to check that the oxygen has at least one hydrogen (i.e., is not esterified).
    for pos in [1, 2]:
        atom_idx = match[pos]
        atom = mol.GetAtomWithIdx(atom_idx)
        # Find the oxygen neighbor that is not the central carbon.
        oxy_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 8]
        if not oxy_neighbors:
            return False, "Terminal CH2 group not attached to oxygen as expected"
        # We expect one such oxygen. (There might be more but we check at least one passes.)
        valid_OH = False
        for o_atom in oxy_neighbors:
            # For a free hydroxyl oxygen, we expect at least one hydrogen.
            # GetTotalNumHs(includeNeighbors=True) returns the number of hydrogens.
            if o_atom.GetTotalNumHs() >= 1:
                valid_OH = True
                break
        if not valid_OH:
            return False, "One or both terminal oxygens are not free hydroxyls; possible extra acylation detected"
    
    return True, "Molecule contains a glycerol backbone with a single acyl ester at the 2-position"

# For testing purposes (can be removed in production):
if __name__ == "__main__":
    test_smiles = [
        "O=C(CCC/C=C\\C\\C=C/C=C/C(C/C=C\\CCCCC)OO)OC(CO)CO",  # 12-HPETE 2-glyceryl ester
        "CCCC\\C=C/C\\C=C/CCCCCCCC(=O)OC(CO)CO",                # 2-linoleoylglycerol
        "CCCCCCCC\\C=C/CCCCCCCC(=O)OC(CO)CO",                  # 2-oleoylglycerol
        "CCCCCCCCCCCCCCCC(=O)OC(CO)CO",                        # 2-palmitoylglycerol
        "[H]C(CO)(CO)OC(=O)CCC\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCCC",  # 2-arachidonoylglycerol
        # A likely false positive (triacylglycerol fragment) should be rejected:
        "O(C(=O)CCCCCCCCCCCCCC)[C@@H](COC(=O)CCCCCCC/C=C\\CCCCCC)COC(=O)CCCC/C=C\\C/C=C\\C/C=C\\C/C=C\\CC",
    ]
    
    for s in test_smiles:
        flag, reason = is_2_monoglyceride(s)
        print(f"SMILES: {s}\nResult: {flag}\nReason: {reason}\n")