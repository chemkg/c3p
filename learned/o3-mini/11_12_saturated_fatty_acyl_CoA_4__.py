"""
Classifies: CHEBI:84948 11,12-saturated fatty acyl-CoA(4-)
"""
#!/usr/bin/env python3
"""
Classifies: CHEBI:84947 – 11,12-saturated fatty acyl-CoA(4-)
Definition:
   Any fatty acyl-CoA(4-) in which the 11-12 bond of the fatty acyl group is saturated.
   (Carbon numbering for the acyl group starts at the carbonyl carbon, which is C1.)
Improved strategy:
  1. Verify the molecule bears a CoA-like moiety by matching both a thioester fragment ([CX3](=O)[S])
     and a nucleotide/adenosine substructure (n1cnc2c(ncnc12)).
  2. Identify the thioester group and consider its carbonyl carbon as the start (C1) of the fatty acyl chain.
  3. From the carbonyl, perform a DFS search – over only carbon atoms – to extract the longest simple path.
  4. Confirm that the chain is long enough (at least 12 carbons) to allow inspection of the 11-12 bond.
  5. Check that the obtained chain is “linear” (i.e. not branched) by ensuring that for each internal carbon
     no additional C–C bond is present beyond the two that form the main chain.
  6. Finally, verify that the bond between C11 and C12 (i.e. indices 10 and 11 in a 0-indexed list)
     is a single (saturated) bond.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_11_12_saturated_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines whether the given SMILES corresponds to an 11,12-saturated fatty acyl-CoA(4-).

    Steps:
       - Parse the molecule.
       - Verify that the molecule contains both a thioester linkage and an adenosine-like (purine) moiety.
       - Identify the thioester match and assume its carbonyl carbon (first atom of the match) is C1.
       - From the thioester, use a depth-first search (DFS) to find the longest continuous carbon chain.
       - Ensure the chain is long enough (>=12 carbons) for an 11-12 bond.
       - Check that the chain is linear (each internal carbon is linked only to its immediate neighbors in the chain).
       - Confirm that the bond between atoms 11 and 12 (chain indices 10 and 11) is a single (saturated) bond.

    Args:
       smiles (str): SMILES string of the molecule.

    Returns:
       bool: True if the molecule qualifies as an 11,12-saturated fatty acyl-CoA(4-), else False.
       str: Reason for the classification.
    """
    # Parse SMILES:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # --- 1. Check for a CoA-like moiety ---
    # a) Look for a thioester bond: [CX3](=O)[S] should capture the acyl-CoA connection.
    thioester_pattern = Chem.MolFromSmarts("[CX3](=O)[S]")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "Thioester (acyl-CoA linkage) not found"

    # b) Look for an adenosine-like purine fragment that is found in CoA:
    nucleotide_pattern = Chem.MolFromSmarts("n1cnc2c(ncnc12)")
    if not mol.HasSubstructMatch(nucleotide_pattern):
        return False, "CoA nucleotide moiety not found"

    # --- 2. Identify the thioester match and treat the carbonyl carbon as C1
    # We assume the first match is the fatty acyl connection.
    thioester_match = thioester_matches[0]
    # In our SMARTS, the first atom is the carbonyl carbon (defined as C1)
    acyl_carbon_idx = thioester_match[0]
    acyl_carbon = mol.GetAtomWithIdx(acyl_carbon_idx)
    
    # Identify neighbor(s) (must be carbon) that lead into the acyl chain.
    acyl_chain_neighbors = [nbr for nbr in acyl_carbon.GetNeighbors() if nbr.GetAtomicNum() == 6]
    if not acyl_chain_neighbors:
        return False, "Fatty acyl chain not found (no carbon attached to acyl carbon)"

    # --- 3. DFS to find the longest simple carbon chain from the carbonyl (C1)
    # DFS function returns a list of atom indices in the path starting at a given atom.
    def dfs(atom, visited):
        best_path = [atom.GetIdx()]
        for nbr in atom.GetNeighbors():
            # Only consider carbon atoms:
            if nbr.GetAtomicNum() != 6:
                continue
            # Do not revisit atoms already in path:
            if nbr.GetIdx() in visited:
                continue
            path = dfs(nbr, visited | {nbr.GetIdx()})
            if len(path) + 1 > len(best_path):
                best_path = [atom.GetIdx()] + path
        return best_path

    longest_chain = []
    # Try from each neighbor and pick the longest resulting chain (prepend the carbonyl atom)
    for neighbor in acyl_chain_neighbors:
        path = dfs(neighbor, {acyl_carbon.GetIdx(), neighbor.GetIdx()})
        full_path = [acyl_carbon.GetIdx()] + path
        if len(full_path) > len(longest_chain):
            longest_chain = full_path

    if len(longest_chain) < 12:
        return False, f"Fatty acyl chain too short (length {len(longest_chain)}); at least 12 carbons required for 11-12 bond check"
    
    # --- 4. Check that the chain is linear (unbranched) ---
    # For the first atom (C1) and the terminal atom, the connectivity in the chain is allowed.
    # For each internal atom, ensure that among its carbon neighbors, only the two from the main chain appear.
    chain_set = set(longest_chain)
    # Get neighbors in the chain for a given atom index (if in the chain)
    for i, idx in enumerate(longest_chain):
        atom = mol.GetAtomWithIdx(idx)
        # Determine expected count: for first and last, expect 1; for internal atoms, exactly 2.
        expected = 1 if i==0 or i==len(longest_chain)-1 else 2
        chain_neighbors = 0
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6 and nbr.GetIdx() in chain_set:
                chain_neighbors += 1
        if chain_neighbors != expected:
            return False, f"Fatty acyl chain appears branched (atom index {idx} has {chain_neighbors} chain-connected carbons, expected {expected})"

    # --- 5. Check that the 11-12 bond (chain positions 10 and 11) is a single bond ---
    idx_C11 = longest_chain[10]
    idx_C12 = longest_chain[11]
    bond_11_12 = mol.GetBondBetweenAtoms(idx_C11, idx_C12)
    if bond_11_12 is None:
        return False, "Bond between C11 and C12 not found"
    if bond_11_12.GetBondType() != Chem.BondType.SINGLE:
        return False, "11-12 bond is not a single (saturated) bond"

    return True, "Molecule is an 11,12-saturated fatty acyl-CoA(4-) with a proper (linear) fatty acyl chain and CoA moiety"


# ---------------------------
# For testing purposes (this block may be removed in production)
if __name__ == "__main__":
    test_smiles = [
        # True positive: 3-oxoisooctadecanoyl-CoA(4-)
        r"CC(C)CCCCCCCCCCCCC(=O)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)n1cnc2c(N)ncnc12",
        # False positive example should be rejected because even if the 11-12 bond is single, the chain is branched.
        r"CC(C)CCCC(C)CCCC(C)CCCC(C)C(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12",
        # Another molecule with a double bond at the 11-12 position
        r"CCCC\C=C/CCCCCCCC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)n1cnc2c(N)ncnc12"
    ]
    for s in test_smiles:
        res, reason = is_11_12_saturated_fatty_acyl_CoA_4__(s)
        print(res, reason)