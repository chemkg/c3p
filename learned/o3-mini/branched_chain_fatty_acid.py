"""
Classifies: CHEBI:35819 branched-chain fatty acid
"""
#!/usr/bin/env python
"""
Classifies: branched-chain fatty acid (BCFA).

A branched-chain fatty acid (BCFA) is a fatty acid (a carboxylic acid with a long aliphatic acyl chain)
that has one or more alkyl substituents (typically a methyl or similarly short branch) on its main chain.
This implementation:
  • Rejects molecules that contain atoms besides carbon and oxygen.
  • Checks that the oxygen count is typical for a fatty acid (either 2 or 3 oxygens).
  • Locates a single carboxyl group (–COOH or deprotonated –COO–).
  • Identifies the alpha carbon (the first carbon attached to the carboxyl carbon).
  • Builds a graph of acyclic (non‐ring) carbon atoms (the fatty acyl chain is assumed acyclic).
  • Uses depth‐first search (DFS) from the alpha carbon to find the longest “main chain”.
  • For each carbon in that main chain, examines neighboring acyclic carbons not on the main chain.
    If such a “branch” is detected and its total number of carbons is short (≤ 3 atoms), it is
    considered a valid alkyl substituent, and the molecule classifies as a BCFA.
    
Note: This is a heuristic approach; different thresholds (minimum main chain length, branch length)
and additional filters (such as rejecting molecules with multiple carboxy groups) can be applied.
"""

from rdkit import Chem

def is_branched_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a branched-chain fatty acid (BCFA) based on its SMILES string.
    
    The algorithm:
      1. Parses the molecule and rejects those with atoms besides C and O.
      2. Verifies that the overall oxygen count is what we expect (2 or 3 atoms).
      3. Uses SMARTS to locate a single carboxyl group (-COOH or its deprotonated form).
      4. Finds the alpha carbon (first non-ring carbon attached to the carboxyl carbon).
      5. Constructs a connectivity graph of all acyclic carbon atoms.
      6. Uses DFS to determine the longest acyclic chain (the “main chain”) starting at the alpha carbon.
      7. For each main-chain carbon, any neighboring acyclic carbons that are not part
         of the main chain are explored. If the branch (its total number of carbon atoms)
         is small (≤ 3 atoms), then it is considered an alkyl substituent, and the molecule
         is declared a BCFA.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule is classified as a branched-chain fatty acid, else False.
        str: A brief explanation for the classification decision.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Step 1: Only allow carbon and oxygen (hydrogens are implicit).
    allowed_atomic_nums = {6, 8}
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in allowed_atomic_nums:
            return False, f"Contains atom {atom.GetSymbol()} (atomic number {atom.GetAtomicNum()}), not typical for a fatty acid"
    
    # Step 2: Check total oxygen count.
    ox_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if ox_count not in {2, 3}:
        return False, f"Oxygen count of {ox_count} not typical for a fatty acid (expected 2 or 3)"
    
    # Step 3: Find the carboxyl group.
    # This SMARTS covers both -COOH and deprotonated -COO–.
    carboxyl_smarts = "[CX3](=O)[OX1H0-,OX2H1]"
    carboxyl_query = Chem.MolFromSmarts(carboxyl_smarts)
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_query)
    if not carboxyl_matches:
        return False, "No carboxyl group found; not a fatty acid"
    if len(carboxyl_matches) != 1:
        return False, "Multiple carboxyl groups detected; not a typical fatty acid"
    # In the SMARTS, the first atom is the carboxyl (carbonyl) carbon.
    carboxyl_C_idx = carboxyl_matches[0][0]
    carboxyl_atom = mol.GetAtomWithIdx(carboxyl_C_idx)
    
    # Step 4: Identify the alpha carbon.
    alpha_carbon_idx = None
    for nbr in carboxyl_atom.GetNeighbors():
        if nbr.GetAtomicNum() == 6 and not nbr.IsInRing():
            alpha_carbon_idx = nbr.GetIdx()
            break
    if alpha_carbon_idx is None:
        return False, "No suitable alpha carbon attached to the carboxyl group; not a fatty acid"
    
    # Step 5: Build the acyclic carbon graph.
    # We only consider carbon atoms that are not in a ring.
    acyclic_carbon_indices = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6 and not atom.IsInRing()]
    if not acyclic_carbon_indices:
        return False, "No acyclic carbon atoms found"
    acyclic_carbon_set = set(acyclic_carbon_indices)
    # Build a dictionary mapping each acyclic carbon index to its neighboring acyclic carbons.
    carbon_graph = { idx: [] for idx in acyclic_carbon_set }
    for idx in acyclic_carbon_set:
        atom = mol.GetAtomWithIdx(idx)
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6 and nbr.GetIdx() in acyclic_carbon_set:
                carbon_graph[idx].append(nbr.GetIdx())
    
    # Step 6: Determine the main (longest) acyclic chain starting at the alpha carbon.
    # We use DFS (disallowing revisiting nodes) to find the longest path.
    def dfs_longest(current, visited):
        longest_path = [current]
        for nbr in carbon_graph.get(current, []):
            if nbr in visited:
                continue
            candidate_path = [current] + dfs_longest(nbr, visited | {nbr})
            if len(candidate_path) > len(longest_path):
                longest_path = candidate_path
        return longest_path

    main_chain = dfs_longest(alpha_carbon_idx, {alpha_carbon_idx})
    if len(main_chain) < 2:
        return False, "Alkyl chain too short to be considered a fatty acid"
    
    # Step 7: For every atom in the main chain, check for branch substituents.
    # A valid branch here is a connected set of acyclic carbons, disjoint from the main chain,
    # with a small total number of carbons (here, we require size <= 3, e.g. methyl, ethyl, or propyl).
    main_chain_set = set(main_chain)
    # Helper: given a starting branch atom, get all connected branch carbons (disallowing main-chain atoms).
    def get_branch_atoms(start):
        branch_atoms = set()
        stack = [start]
        while stack:
            cur = stack.pop()
            if cur in branch_atoms:
                continue
            branch_atoms.add(cur)
            for nb in carbon_graph.get(cur, []):
                if nb in main_chain_set:
                    continue
                if nb not in branch_atoms:
                    stack.append(nb)
        return branch_atoms

    for idx in main_chain:
        for nbr in carbon_graph.get(idx, []):
            if nbr not in main_chain_set:
                branch_atoms = get_branch_atoms(nbr)
                branch_size = len(branch_atoms)
                # For a typical BCFA, the branch (alkyl substituent) is very short.
                if branch_size <= 3:
                    return True, f"Branching detected at main chain carbon atom index {idx} (branch size: {branch_size})"
    
    return False, "No suitable short branching substituents detected on the fatty acyl chain"

# Example usage (for debugging/testing):
if __name__ == '__main__':
    test_smiles = [
        "OC(=O)CC(C)=C",                         # isopropenylacetic acid: expected BCFA
        "CCC(C)CCCCCCCCCCCCCCCCCCCCCCCCCCC(O)=O", # 28-methyltriacontanoic acid: expected BCFA
        "OC(=O)CCCCCCCC",                         # linear fatty acid: expected NOT BCFA
        "OC(=O)C(CC)CC",                          # 2-Ethylbutanoic acid: expected BCFA
    ]
    for s in test_smiles:
        result, reason = is_branched_chain_fatty_acid(s)
        print(s, "->", result, ":", reason)