"""
Classifies: CHEBI:35819 branched-chain fatty acid
"""
"""
Classifies: branched-chain fatty acid.
A branched-chain fatty acid (BCFA) is a fatty acid (a carboxylic acid with an aliphatic chain)
in which one or more alkyl groups (typically simple carbon‐only substituents) branch off from its main chain.
This implementation:
  1. Checks for a carboxyl group (–COOH or –COO–).
  2. Locates the alpha carbon (the first carbon attached to the carboxyl carbon).
  3. Constructs a “main” acyl chain by a depth‐first search (DFS) on carbon atoms—restricting to acyclic (non‐ring) carbons.
  4. Examines each atom in the main chain for a carbon neighbor (apart from the expected connections)
     and then checks that that branch is “simple” (i.e. comprised solely of carbon atoms).
If found, the molecule is classified as a branched‐chain fatty acid.
"""

from rdkit import Chem

def is_branched_chain_fatty_acid(smiles: str):
    """
    Determines whether a molecule is a branched-chain fatty acid (BCFA) from its SMILES string.
    
    The algorithm:
      1. Parses the molecule and checks for a carboxyl moiety using a SMARTS.
      2. Finds the alpha carbon attached to the carboxyl carbon.
      3. Constructs the longest acyclic carbon chain (excluding ring atoms) from the alpha carbon.
      4. For each carbon in this main chain, any extra (carbon) neighbor not in the chain is
         examined to ensure it is a “simple” alkyl branch (i.e. contains only carbons).
      5. If any such valid branch is found, the molecule is declared a BCFA.
    
    Args:
      smiles (str): SMILES string for the molecule.
      
    Returns:
      (bool, str): Tuple with classification (True if BCFA) and a reason message.
    """
    # Parse molecule from SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # SMARTS for a carboxyl group: matches either -COOH or deprotonated -COO-
    carboxyl_smarts = "[CX3](=O)[OX1H0-,OX2H1]"
    carboxyl_query = Chem.MolFromSmarts(carboxyl_smarts)
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_query)
    if not carboxyl_matches:
        return False, "No carboxyl group found; not a fatty acid"
    
    # Take the first carboxyl match. In our SMARTS the first atom is the carboxyl carbon.
    carboxyl_atoms = carboxyl_matches[0]
    carboxyl_C = carboxyl_atoms[0]
    
    # Find the alpha carbon: a carbon neighbor of the carboxyl carbon.
    carboxyl_atom = mol.GetAtomWithIdx(carboxyl_C)
    alpha_carbon = None
    for nbr in carboxyl_atom.GetNeighbors():
        if nbr.GetAtomicNum() == 6:  # must be a carbon
            alpha_carbon = nbr.GetIdx()
            break
    if alpha_carbon is None:
        return False, "No alpha carbon attached to carboxyl group; not a fatty acid"
    
    # Build a set of all carbon atom indices in the molecule.
    carbon_indices = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    # Also build a mapping: atom index -> list of neighboring carbon indices.
    carbon_neighbors = {idx: [] for idx in carbon_indices}
    for idx in carbon_indices:
        atom = mol.GetAtomWithIdx(idx)
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6:
                carbon_neighbors[idx].append(nbr.GetIdx())
    
    # A helper DFS function that finds the longest chain (as a list of atom indices)
    # starting from the given index. We restrict the search to carbon atoms that are not in rings.
    def dfs_longest(current, path):
        longest = path
        for nbr in carbon_neighbors.get(current, []):
            # Skip if already in path or if the neighbor is in a ring.
            if nbr in path:
                continue
            if mol.GetAtomWithIdx(nbr).IsInRing():
                continue
            new_path = dfs_longest(nbr, path + [nbr])
            if len(new_path) > len(longest):
                longest = new_path
        return longest

    # Get the main chain starting at the alpha carbon.
    main_chain = dfs_longest(alpha_carbon, [alpha_carbon])
    if len(main_chain) < 1:
        return False, "Alkyl chain too short to be a fatty acid"
    
    # For additional filtering: if any atom in the main chain is in a ring, reject.
    for idx in main_chain:
        if mol.GetAtomWithIdx(idx).IsInRing():
            return False, "Main chain contains ring atoms; not a typical fatty acyl chain"
    
    # Define a helper function to check if a branch (a connected component starting at branch_idx)
    # is a "simple alkyl substituent" (i.e. all atoms accessible are carbon) while not traversing
    # into the main acyl chain.
    def is_simple_alkyl(branch_idx, exclude_set):
        visited = set()
        stack = [branch_idx]
        while stack:
            cur = stack.pop()
            if cur in visited:
                continue
            visited.add(cur)
            atom = mol.GetAtomWithIdx(cur)
            # Check every neighbor (full connectivity)
            for nbr in atom.GetNeighbors():
                nbr_idx = nbr.GetIdx()
                # Skip neighbors that are part of the main chain (or any explicitly excluded atom)
                if nbr_idx in exclude_set:
                    continue
                # If any neighbor is not carbon, then this branch is not purely alkyl.
                if nbr.GetAtomicNum() != 6:
                    return False
                if nbr_idx not in visited:
                    stack.append(nbr_idx)
        return True

    # Mark the main chain as a set for quick lookup.
    main_chain_set = set(main_chain)
    
    # Examine each carbon in the main chain for extra carbon substituents.
    # For the first atom in the chain (the alpha carbon), allow the connection back to the carboxyl carbon.
    for idx in main_chain:
        atom = mol.GetAtomWithIdx(idx)
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() != 6:
                continue
            nbr_idx = nbr.GetIdx()
            # For the alpha carbon, allow its bond back to the carboxyl carbon.
            if idx == main_chain[0] and nbr_idx == carboxyl_C:
                continue
            # If this neighbor is not in the main chain it is a candidate branch.
            if nbr_idx not in main_chain_set:
                # Check if the branch is a simple alkyl group.
                if is_simple_alkyl(nbr_idx, main_chain_set):
                    return True, f"Branching detected at carbon atom index {idx}"
    
    return False, "No branching substituents detected on the fatty acyl chain"

# Example usage (uncomment to test):
# smiles_list = [
#     "OC(=O)CC(C)=C",  # isopropenylacetic acid, expected BCFA
#     "OC(=O)C(CC)CC",  # 2-Ethylbutanoic acid, expected BCFA
#     "CCC(C)CCCCCCCCCCCCCCCCCCCCCCCCCCC(O)=O",  # 28-methyltriacontanoic acid, expected BCFA
# ]
# for s in smiles_list:
#     result, reason = is_branched_chain_fatty_acid(s)
#     print(s, "->", result, reason)