"""
Classifies: CHEBI:25413 monounsaturated fatty acid
"""
"""
Classifies: Monounsaturated fatty acid.
Definition: A fatty acid (free acid or esterified acyl chain) that contains exactly one
carbon–carbon unsaturation (double or triple bond) in the acyl chain. In our implementation,
we look for a carboxyl/ester group and then for each candidate carbon chain neighbor (non-oxygen)
we enumerate all linear acyl chain paths (only aliphatic, non-aromatic, non‐cyclic carbons). 
We then select the longest path and verify that it contains exactly one unsaturated bond.
MUFAs are known to have positive effects on the cardiovascular system and in treatment of diabetes.
"""

from rdkit import Chem

def is_monounsaturated_fatty_acid(smiles: str):
    """
    Determines if a molecule contains a monounsaturated fatty acyl chain fragment.
    For each acidic (or ester) carbonyl group, every eligible carbon neighbor is considered as
    a potential starting point of an acyl chain. Then, using a DFS that only follows aliphatic,
    non-aromatic, non-cyclic carbon atoms, every simple path (i.e. linear fragment) is enumerated.
    The longest path is chosen as representing the acyl chain. If that path has length >= 2 (in 
    terms of number of carbons) and contains exactly one carbon–carbon unsaturation (double or triple bond),
    then the molecule is classified as containing a monounsaturated fatty acid.
    
    Args:
      smiles (str): SMILES string representing the molecule.
      
    Returns:
      bool: True if molecule contains at least one fatty acyl chain meeting the MUFA criteria.
      str: Explanation for the classification decision.
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS pattern for free acid (-COOH) and ester (-COO-)
    free_acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    # The ester pattern here requires the oxygen to be attached to a carbon 
    # (it may be embedded in a larger structure such as glycerol or phospholipid)
    ester_pattern = Chem.MolFromSmarts("[CX3](=O)O[#6]")
    
    free_acid_matches = mol.GetSubstructMatches(free_acid_pattern)
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    
    # Allow both types of matches
    candidate_matches = list(free_acid_matches) + list(ester_matches)
    if not candidate_matches:
        return False, "No free acid or ester carbonyl group found to derive a fatty acyl chain."
    
    # For each match, we try every carbon neighbor of the carbonyl (skipping the paired oxygen)
    # as the potential start of a fatty acyl chain.
    candidates = []  # each candidate is (carbonyl_idx, start_idx)
    for match in candidate_matches:
        carbonyl_idx = match[0]  # assumed to be the carbonyl carbon
        # Get all neighbors of the carbonyl atom that are carbon (atomic number 6).
        carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
        for nbr in carbonyl_atom.GetNeighbors():
            # Exclude oxygen neighbors (which are part of -COOH or -COO-)
            if nbr.GetAtomicNum() == 6:
                candidates.append((carbonyl_idx, nbr.GetIdx()))
    
    if not candidates:
        return False, "No valid candidate fatty acyl chain found (carbonyl center not attached to any carbon)."

    # Define a helper function to enumerate all simple (linear) paths from a given starting atom.
    # We only allow traversal through aliphatic carbons (atomic number 6) that are neither aromatic nor in a ring.
    def dfs_paths(current_path):
        last_idx = current_path[-1]
        last_atom = mol.GetAtomWithIdx(last_idx)
        extended = False
        for nbr in last_atom.GetNeighbors():
            nbr_idx = nbr.GetIdx()
            # continue only if this neighbor is an aliphatic carbon and not already in path
            if nbr_idx in current_path:
                continue
            if nbr.GetAtomicNum() != 6 or nbr.GetIsAromatic() or nbr.IsInRing():
                continue
            # Extend the path.
            extended = True
            yield from dfs_paths(current_path + [nbr_idx])
        # If we could not extend further, yield the current path as maximal.
        if not extended:
            yield current_path

    # For each candidate, enumerate all possible acyl chain linear paths.
    for carbonyl_idx, start_idx in candidates:
        # Only consider the starting atom if it is an eligible carbon.
        start_atom = mol.GetAtomWithIdx(start_idx)
        if start_atom.GetAtomicNum() != 6 or start_atom.GetIsAromatic() or start_atom.IsInRing():
            continue
        all_paths = list(dfs_paths([start_idx]))
        if not all_paths:
            continue
        # Choose the longest linear path (in terms of number of carbon atoms in the chain)
        longest_path = max(all_paths, key=lambda p: len(p))
        if len(longest_path) < 2:
            continue  # chain too short; try next candidate
        
        # Count the number of unsaturated carbon–carbon bonds along the selected path.
        unsat_count = 0
        # Note: we only consider bonds between consecutive atoms in the path.
        for i in range(len(longest_path) - 1):
            bond = mol.GetBondBetweenAtoms(longest_path[i], longest_path[i+1])
            if bond is None:
                continue
            if bond.GetBondType() in (Chem.rdchem.BondType.DOUBLE, Chem.rdchem.BondType.TRIPLE):
                unsat_count += 1
        
        # We require that the acyl chain has exactly one carbon–carbon unsaturation.
        if unsat_count == 1:
            return True, ("Found acyl chain (with %d carbon atoms) derived from a carbonyl at atom %d having exactly one "
                          "carbon–carbon unsaturation." % (len(longest_path), carbonyl_idx))
    # If none of the candidate chains qualify, return False with details.
    return False, ("No fatty acyl chain candidate with exactly one unsaturation was found among %d candidate(s)." 
                   % (len(candidates)))

# Example usage (uncomment to test a few examples):
# examples = [
#     ("CC(C)=CC(O)=O", "3-methylbut-2-enoic acid"),
#     ("OCCCCCCCCCCCCCCC\\C=C\\C(O)=O", "(2E)-18-hydroxyoctadec-2-enoic acid"),
#     ("CCCCCC\\C=C\\CCCCCCCCCC(O)=O", "trans-vaccenic acid"),
#     ("OC(=O)CCCC/C=C/CC", "(2E)-oct-2-enoic acid")
# ]
# for smi, name in examples:
#     result, reason = is_monounsaturated_fatty_acid(smi)
#     print(name, "->", result, reason)