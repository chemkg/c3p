"""
Classifies: CHEBI:35819 branched-chain fatty acid
"""
#!/usr/bin/env python
"""
Classifies: branched-chain fatty acid (BCFA).

A branched-chain fatty acid is a fatty acid (a carboxylic acid with an aliphatic acyl chain)
which has one or more alkyl substituents (typically a methyl branch or similar) on its main chain.
This implementation:
  1. Rejects molecules containing atoms other than carbon and oxygen.
  2. Detects a carboxyl group (–COOH or deprotonated –COO–) via SMARTS.
  3. Finds the alpha carbon (the acyclic carbon directly attached to the carboxyl carbon).
  4. Restricts the analysis to acyclic (non‐ring) carbon atoms (the “fatty acyl chain”).
  5. Finds the longest acyclic chain (by DFS) starting from the alpha carbon.
  6. Checks each main‐chain carbon for an extra (branch) neighbor (a carbon not in the main chain).
If any valid branch is found, the molecule is declared a branched‐chain fatty acid.
"""

from rdkit import Chem

def is_branched_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a branched-chain fatty acid (BCFA) based on its SMILES string.
    
    The algorithm:
      1. Parses the molecule and rejects those with heteroatoms (atoms besides C and O).
      2. Uses SMARTS to locate a carboxyl group (-COOH or -COO–).
      3. Finds the alpha carbon attached to the carboxyl carbon.
      4. Constructs an explicit connectivity graph of acyclic (non‐ring) carbon atoms.
      5. Uses a DFS to find the longest acyclic main chain starting at the alpha carbon.
      6. For each carbon in the main chain, any extra carbon not in the main chain is marked as a branch.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule is classified as a branched-chain fatty acid, else False.
        str: Reason for classification.
    """
    # Step 1: Parse molecule and check for allowed elements.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    # For our purpose, we restrict to fatty acids made solely of carbon and oxygen.
    # (Hydrogens are implicit.) This helps reject peptides and other compounds.
    allowed_atomic_nums = {6, 8}
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in allowed_atomic_nums:
            return False, f"Contains atom {atom.GetSymbol()} (atomic number {atom.GetAtomicNum()}), not typical for a fatty acid"
    
    # Step 2: Find a carboxyl group.
    # The following SMARTS covers both -COOH and deprotonated -COO–:
    carboxyl_smarts = "[CX3](=O)[OX1H0-,OX2H1]"
    carboxyl_query = Chem.MolFromSmarts(carboxyl_smarts)
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_query)
    if not carboxyl_matches:
        return False, "No carboxyl group found; not a fatty acid"
    # Use the first matching carboxyl group. In the SMARTS, the first atom is the carbonyl carbon.
    carboxyl_C = carboxyl_matches[0][0]
    carboxyl_atom = mol.GetAtomWithIdx(carboxyl_C)
    
    # Step 3: Find the alpha carbon (first carbon attached to the carboxyl carbon)
    alpha_carbon = None
    for nbr in carboxyl_atom.GetNeighbors():
        if nbr.GetAtomicNum() == 6 and not nbr.IsInRing():
            alpha_carbon = nbr.GetIdx()
            break
    if alpha_carbon is None:
        return False, "No suitable alpha carbon attached to the carboxyl group; not a fatty acid"
    
    # Step 4: Build a graph of acyclic carbons in the molecule.
    # We consider only carbon atoms that are not in any ring.
    carbon_indices = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6 and not atom.IsInRing()]
    if not carbon_indices:
        return False, "No acyclic carbon atoms found"
    acyclic_carbon_set = set(carbon_indices)
    # Construct a dictionary mapping each acyclic carbon index to its neighboring acyclic carbons.
    carbon_graph = { idx: [] for idx in acyclic_carbon_set }
    for idx in acyclic_carbon_set:
        atom = mol.GetAtomWithIdx(idx)
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6 and nbr.GetIdx() in acyclic_carbon_set:
                carbon_graph[idx].append(nbr.GetIdx())
    
    # Step 5: Starting at the alpha carbon, perform DFS to find the longest linear chain.
    # We disallow revisiting nodes.
    def dfs_longest(current, visited):
        longest_path = [current]
        for nbr in carbon_graph.get(current, []):
            if nbr in visited:
                continue
            path = dfs_longest(nbr, visited | {nbr})
            # Prepend current to the candidate path
            candidate = [current] + path
            if len(candidate) > len(longest_path):
                longest_path = candidate
        return longest_path

    main_chain = dfs_longest(alpha_carbon, {alpha_carbon})
    # A minimal fatty acid should have at least 2 carbons in its chain (alpha + one more)
    if len(main_chain) < 2:
        return False, "Alkyl chain too short to be considered a fatty acid"
    
    # Step 6: Check each carbon in the main chain for extra branch substituents.
    main_chain_set = set(main_chain)
    # For each atom in the main chain, see if it has a neighboring acyclic carbon that is not part of the main chain.
    for idx in main_chain:
        for nbr in carbon_graph.get(idx, []):
            if nbr not in main_chain_set:
                # We found a candidate branch. Perform a simple DFS on the branch (excluding the main chain)
                branch_atoms = set()
                stack = [nbr]
                while stack:
                    current_branch = stack.pop()
                    if current_branch in branch_atoms:
                        continue
                    branch_atoms.add(current_branch)
                    for nb in carbon_graph.get(current_branch, []):
                        if nb in main_chain_set:
                            continue
                        if nb not in branch_atoms:
                            stack.append(nb)
                if branch_atoms:
                    return True, f"Branching detected at main chain carbon atom index {idx}"
    
    return False, "No branching substituents detected on the fatty acyl chain"

# Example usage (uncomment to test):
# test_smiles = [
#     "OC(=O)CC(C)=C",                        # isopropenylacetic acid: expected BCFA (branch off main chain)
#     "CCC(C)CCCCCCCCCCCCCCCCCCCCCCCCCCC(O)=O",# 28-methyltriacontanoic acid: expected BCFA
#     "OC(=O)CCCCCCCC",                        # linear fatty acid: expected not BCFA
#     "OC(=O)C(CC)CC",                         # 2-Ethylbutanoic acid: expected BCFA
# ]
# for s in test_smiles:
#     result, reason = is_branched_chain_fatty_acid(s)
#     print(s, "->", result, ":", reason)