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
  1. Check for CoA by requiring both a thioester linkage ([CX3](=O)[S]) and a purine-like fragment 
     (n1cnc2c(ncnc12)) that is present in CoA(4-).
  2. Starting from the thioester (carbonyl carbon), use a DFS search to find the longest linear chain 
     of carbons (thus capturing the main fatty acyl chain even with small branches).
  3. Confirm that the chain is long enough (>= 12 carbons) and that the bond between carbon 11 and 
     carbon 12 is a single (saturated) bond.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_11_12_saturated_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines whether the given SMILES string corresponds to an 11,12-saturated fatty acyl-CoA(4-).
    
    Steps:
       - Check that the molecule is valid.
       - Verify that it contains a CoA-like fragment (using a thioester SMARTS and a purine-ring SMARTS).
       - Locate the thioester group. Assume the carbonyl carbon in the match is C1.
       - From the thioester, perform a DFS search (over carbon atoms only) to find the longest simple carbon chain.
       - Verify that the chain has at least 12 carbons and check that specifically the bond between 
         carbon 11 and 12 (i.e. chain indices 10 and 11 with C1 as index 0) is a single bond.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule qualifies as an 11,12-saturated fatty acyl-CoA(4-), False otherwise.
        str: Explanation message.
    """
    # Parse the SMILES:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # ---------------------------
    # 1. Check for CoA-like fragment.
    # We require both:
    #   a) A thioester bond: carbonyl attached to sulfur.
    thioester_pattern = Chem.MolFromSmarts("[CX3](=O)[S]")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "Thioester (acyl-CoA linkage) not found"
    #   b) A purine-like ring that is part of CoA's adenosine moiety.
    nucleotide_pattern = Chem.MolFromSmarts("n1cnc2c(ncnc12)")
    if not mol.HasSubstructMatch(nucleotide_pattern):
        return False, "CoA nucleotide moiety not found"
    
    # ---------------------------
    # 2. Identify the thioester match and choose the first match.
    # We assume that the first match corresponds to the fatty acyl-CoA connection.
    thioester_match = thioester_matches[0]
    # In our SMARTS, the first atom is the carbonyl carbon (C1, the start of the fatty acyl chain)
    acyl_carbon_idx = thioester_match[0]
    acyl_carbon = mol.GetAtomWithIdx(acyl_carbon_idx)
    
    # Find the neighbor in the acyl chain: choose among neighbors the one that is carbon (atomic num 6)
    acyl_chain_neighbors = [nbr for nbr in acyl_carbon.GetNeighbors() if nbr.GetAtomicNum() == 6]
    if not acyl_chain_neighbors:
        return False, "Fatty acyl chain not found (no carbon attached to acyl carbon)"
    # For our purposes, if more than one neighbor, we will search all starting paths.
    
    # ---------------------------
    # 3. Find the longest carbon chain (a simple path) starting from the acyl chain.
    # Use DFS that collects simple paths (no repeated atoms) and, among those starting at the neighbor,
    # we choose the longest path. Then we prepend the thioester carbon.
    def dfs(atom, visited):
        """Recursive DFS returning the longest simple path (list of atom indices) starting at atom."""
        best_path = [atom.GetIdx()]
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() != 6:
                continue
            if nbr.GetIdx() in visited:
                continue
            path = dfs(nbr, visited | {nbr.GetIdx()})
            # Prepend current neighbor to path.
            if len(path) + 1 > len(best_path):
                best_path = [atom.GetIdx()] + path
        return best_path

    longest_chain = []
    for neighbor in acyl_chain_neighbors:
        path = dfs(neighbor, {acyl_carbon.GetIdx(), neighbor.GetIdx()})
        # Include the starting acyl carbon (C1) at the beginning.
        full_path = [acyl_carbon.GetIdx()] + path
        if len(full_path) > len(longest_chain):
            longest_chain = full_path

    if len(longest_chain) < 12:
        return False, f"Fatty acyl chain too short (length {len(longest_chain)}); at least 12 carbons required for 11-12 bond check"
    
    # ---------------------------
    # 4. Check that the 11-12 bond is saturated.
    # We number carbons such that the carbonyl carbon is C1, so in 0-indexing, the bond between 
    # chain[10] and chain[11] is the 11-12 bond.
    idx_C11 = longest_chain[10]
    idx_C12 = longest_chain[11]
    bond_11_12 = mol.GetBondBetweenAtoms(idx_C11, idx_C12)
    if bond_11_12 is None:
        return False, "Bond between C11 and C12 not found"
    if bond_11_12.GetBondType() != Chem.BondType.SINGLE:
        return False, "11-12 bond is not a single (saturated) bond"
    
    return True, "Molecule is an 11,12-saturated fatty acyl-CoA(4-) with a proper fatty acyl chain and CoA moiety"


# ---------------------------
# For testing (this block can be removed in production)
if __name__ == "__main__":
    # Test examples:
    test_smiles = [
        # True positive: 3-oxoisooctadecanoyl-CoA(4-)
        r"CC(C)CCCCCCCCCCCCC(=O)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)n1cnc2c(N)ncnc12",
        # False positive example should be rejected (e.g. one with a double bond at C11-C12).
        r"CCCC\C=C/CCCCCCCC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)n1cnc2c(N)ncnc12"
    ]
    for s in test_smiles:
        result, reason = is_11_12_saturated_fatty_acyl_CoA_4__(s)
        print(result, reason)