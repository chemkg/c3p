"""
Classifies: CHEBI:15904 long-chain fatty acid
"""
"""
Classifies: Long-chain fatty acid (chain length ranging from C13 to C22)

A long-chain fatty acid is defined as a free carboxylic acid where a terminal –C(=O)OH
(or its deprotonated form) is attached to exactly one carbon chain. The length of the
longest simple carbon path (starting from the carbon attached to the free acid group and
including the acid carbon) must be between 13 and 22 carbons. In contrast to the previous
approach, rings and small branches (e.g. cyclopropane motifs) are tolerated; the longest
simple carbon path is computed via DFS.
"""
from rdkit import Chem

def is_long_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule qualifies as a long-chain fatty acid (chain between 13 and 22 carbons)
    based on its SMILES string.

    The approach is:
      1. Find a free carboxylic acid group by matching a SMARTS pattern that covers both the protonated and deprotonated forms.
      2. Ensure that the candidate acid carbon is attached to exactly one other carbon (the chain origin).
      3. Do a depth-first search on the subgraph made of carbon atoms only (starting at the chain neighbor)
         to find the longest simple carbon path that can be reached. (Here “simple” means no atom is revisited.)
      4. Count the acid carbon in addition to the chain path. If the total number falls between 13 and 22, return True.
         Otherwise, return False with an appropriate reason.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule qualifies as a long-chain fatty acid, False otherwise.
        str: Reason describing the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS to detect a free carboxylic acid group:
    # Matches a sp2 carbon double-bonded to one oxygen and singly bonded to an oxygen that either has an explicit H
    # or a negative charge (thus matching both –COOH and –COO- forms).
    acid_smarts = Chem.MolFromSmarts("[CX3](=O)[O;H1,-1]")
    acid_matches = mol.GetSubstructMatches(acid_smarts)
    if not acid_matches:
        return False, "No free carboxyl group (-C(=O)O or -C(=O)[O-]) found"

    # We may have several acid groups; try each candidate.
    # We require that the acid carbon (the first atom in the SMARTS match) is attached to exactly one carbon (its chain neighbor)
    # Note: the match returns indices in the order of the pattern, so index 0 is the acid carbon.
    for match in acid_matches:
        acid_carbon_idx = match[0]
        acid_carbon = mol.GetAtomWithIdx(acid_carbon_idx)

        # Identify attached carbon atoms (excluding the two oxygens that form the carboxyl group).
        chain_neighbors = []
        for nbr in acid_carbon.GetNeighbors():
            if nbr.GetAtomicNum() == 6:
                chain_neighbors.append(nbr)
        if len(chain_neighbors) != 1:
            continue  # this candidate is not a terminal free acid (e.g. it might be an ester group if >1 carbon)
        chain_start = chain_neighbors[0]

        # Define a recursive DFS function that traverses only carbon atoms.
        def dfs(atom, visited):
            max_length = 1  # count the current atom
            for nbr in atom.GetNeighbors():
                if nbr.GetAtomicNum() != 6:
                    continue  # only consider carbons
                if nbr.GetIdx() in visited:
                    continue
                # Explore this neighbor further
                new_visited = visited | {nbr.GetIdx()}
                path_length = 1 + dfs(nbr, new_visited)
                if path_length > max_length:
                    max_length = path_length
            return max_length

        # Start DFS from chain_start, and initialize visited with the acid carbon to prevent going backward.
        visited = {acid_carbon_idx, chain_start.GetIdx()}
        chain_path_length = dfs(chain_start, visited)
        total_chain_length = 1 + chain_path_length  # include the acid carbon

        # Check if total chain length (number of carbon atoms) is within [13,22].
        if total_chain_length < 13:
            # Provide reasoning showing the chain is too short.
            return False, f"Chain length too short: {total_chain_length} carbons (< 13 required)"
        if total_chain_length > 22:
            return False, f"Chain length too long: {total_chain_length} carbons (> 22 allowed)"
        
        # If we reached here, this candidate acid group has a carbon chain of acceptable length.
        return True, f"Terminal carboxyl group found with longest carbon chain of {total_chain_length} carbons"

    # If none of the acid groups satisfied the criterion.
    return False, "No terminal free acid group found with exactly one attached carbon"

# (Optional) Example test run:
# test_smiles = [
#     ("O(O)[C@H](CCCCC)\\C=C\\CCCCCCCCCC(O)=O", "13R-HpOME(11E)"),
#     ("C([C@@H](O)CCCCCCCCCCCCC\\C=C\\C(O)=O)", "Example fatty acid"),
# ]
# for s, name in test_smiles:
#     result, reason = is_long_chain_fatty_acid(s)
#     print(f"SMILES: {s}\nNAME: {name}\nResult: {result}, Reason: {reason}\n")