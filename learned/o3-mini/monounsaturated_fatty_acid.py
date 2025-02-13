"""
Classifies: CHEBI:25413 monounsaturated fatty acid
"""
"""
Classifies: Monounsaturated fatty acid.
Definition: Any fatty acid (free acid or esterified acyl chain) with exactly one carbon–carbon unsaturation
(double or triple bond) in the fatty acyl chain. In our approach the candidate acyl chain has to be strictly aliphatic
(no aromatic atoms or ring membership) and must derive from a carboxyl (or ester) carbonyl that is attached
to only one carbon chain.
MUFAs are known to have positive effects on the cardiovascular system and in diabetes treatment.
"""

from rdkit import Chem

def is_monounsaturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a monounsaturated fatty acid (MUFA) based on its SMILES string.
    Unlike our previous version, we now also try to find fatty acyl chain fragments embedded in larger molecules,
    such as in phospholipids. For a candidate acyl chain we require:
      1) A carboxyl-like carbonyl center that is either a free acid (-COOH) or an ester (-COO-).
      2) The carbonyl carbon must be attached to exactly one carbon atom on the acyl chain side.
      3) Starting from that carbon, a contiguous (DFS) aliphatic chain is built – we restrict to non-aromatic,
         non-cyclic carbon atoms.
      4) The chain must have a minimal length (at least 2 carbon atoms) and contain exactly one carbon–carbon 
         unsaturation (double or triple bond) aside from the carbonyl bond.
    
    Args:
      smiles (str): SMILES string representing the molecule.
      
    Returns:
      bool: True if molecule contains at least one fatty acyl chain meeting the MUFA criteria.
      str: Explanation for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS pattern for a free acid group: -COOH
    free_acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    # Define SMARTS pattern for an ester group: -COO-
    ester_pattern = Chem.MolFromSmarts("[CX3](=O)O[#6]")
    
    # Get matches from both patterns. Each match returns a tuple:
    # For both patterns, we assume that index 0 is the carbonyl carbon and index 1 is the O (either -OH or in -OR).
    free_acid_matches = mol.GetSubstructMatches(free_acid_pattern)
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    candidate_matches = list(free_acid_matches) + list(ester_matches)
    
    if not candidate_matches:
        return False, "No free acid or ester carbonyl group found to derive a fatty acyl chain."
    
    # List to hold candidate acyl chain starting points.
    # Each candidate is a tuple: (carbonyl_idx, candidate_start_idx)
    candidates = []
    for match in candidate_matches:
        carbonyl_idx = match[0]
        oxy_idx = match[1]  # the oxygen in the carboxyl or ester group
        carbonyl = mol.GetAtomWithIdx(carbonyl_idx)
        # Look for a neighbor that is a carbon (the acyl chain) excluding the matched oxygen.
        neighbor_c_idxs = [nbr.GetIdx() for nbr in carbonyl.GetNeighbors() 
                             if nbr.GetAtomicNum() == 6 and nbr.GetIdx() != oxy_idx]
        if len(neighbor_c_idxs) == 1:
            # We have a single carbon neighbor; use it as the start of our candidate fatty acyl chain.
            candidates.append((carbonyl_idx, neighbor_c_idxs[0]))
        # If there are 0 or more than one carbon neighbors, we skip this candidate as ambiguous.
    
    if not candidates:
        return False, "No valid fatty acyl chain candidate found (carbonyl center not attached to exactly one carbon)."
    
    # For each candidate, perform a DFS to extract a contiguous acyl chain.
    # We will only travel through carbon atoms that are not aromatic and not in a ring.
    for carbonyl_idx, start_idx in candidates:
        chain_atoms = set()
        def dfs(atom_idx, came_from):
            if atom_idx in chain_atoms:
                return
            atom = mol.GetAtomWithIdx(atom_idx)
            # Only allow aliphatic carbons (atomic number 6), not aromatic and not in a ring.
            if atom.GetAtomicNum() != 6 or atom.GetIsAromatic() or atom.IsInRing():
                return
            chain_atoms.add(atom_idx)
            for nbr in atom.GetNeighbors():
                if nbr.GetAtomicNum() == 6 and nbr.GetIdx() != came_from:
                    dfs(nbr.GetIdx(), atom_idx)
        dfs(start_idx, carbonyl_idx)
        
        # Check that the candidate chain is not trivially short.
        if len(chain_atoms) < 2:
            continue  # try the next candidate
        
        # Now, count the carbon–carbon unsaturated bonds (double or triple) within the extracted chain.
        unsaturation_count = 0
        for bond in mol.GetBonds():
            i = bond.GetBeginAtom().GetIdx()
            j = bond.GetEndAtom().GetIdx()
            # Count only bonds where both endpoints are in the chain.
            if i in chain_atoms and j in chain_atoms:
                if bond.GetBondType() in (Chem.rdchem.BondType.DOUBLE, Chem.rdchem.BondType.TRIPLE):
                    unsaturation_count += 1
        
        # We expect exactly one unsaturation in the acyl chain.
        if unsaturation_count == 1:
            return True, ("Found acyl chain (with %d carbon atoms) derived from a carbonyl at atom %d having exactly one "
                          "carbon–carbon unsaturation." % (len(chain_atoms), carbonyl_idx))
    # If none of the candidate chains qualify, return False with details.
    return False, ("No fatty acyl chain candidate with exactly one unsaturation was found among %d candidate(s)." 
                   % len(candidates))

# Example tests: (uncomment to test examples)
# examples = [
#     ("CC(C)=CC(O)=O", "3-methylbut-2-enoic acid"),
#     ("OCCCCCCCCCCCCCCC\\C=C\\C(O)=O", "(2E)-18-hydroxyoctadec-2-enoic acid"),
#     ("CCCCCC\\C=C\\CCCCCCCCCC(O)=O", "trans-vaccenic acid"),
#     ("OC(=O)[C@]1([C@@]2([C@]([C@]3([C@@]4(CC2)C[C@H](C([C@@H]4O)=C)CC3)[H])(CCC1)C)[H])C", "grandiflorolic acid"),
#     ("OC(=O)CCCC/C=C/CC", "(2E)-oct-2-enoic acid")
# ]
# for smi, name in examples:
#     result, reason = is_monounsaturated_fatty_acid(smi)
#     print(name, "->", result, reason)