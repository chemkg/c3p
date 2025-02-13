"""
Classifies: CHEBI:29348 fatty amide
"""
"""
Classifies: Fatty Amide
Definition:
  A monocarboxylic acid amide derived from a fatty acid.
  In our approach, a fatty amide is recognized by an amide bond C(=O)N
  wherein the acyl (C=O) side is “fatty” (i.e. has an aliphatic chain of
  modest length – even allowing branched chains, here a minimum linear-chain length of 4 is required)
  and the substituents on the amide nitrogen are comparatively short 
  (none having a linear aliphatic chain longer than 4).
  
This heuristic was developed in response to many false positives (e.g. complex peptides)
and false negatives (e.g. dolichotheline) seen previously.
"""

from rdkit import Chem

def longest_aliphatic_chain_length(atom, visited):
    """
    Recursively computes the longest linear chain length starting from the given atom.
    We only traverse carbon atoms that are:
      - non-aromatic,
      - not in any ring.
    'visited' is a set of atom indices already seen.
    (We do not require the chain to be perfectly linear – we take the longest branch.)
    Returns an integer count (at least 1 for the starting atom).
    """
    max_len = 1  # count the current atom
    for nbr in atom.GetNeighbors():
        if nbr.GetAtomicNum() == 6 and (not nbr.GetIsAromatic()) and (not nbr.IsInRing()):
            if nbr.GetIdx() in visited:
                continue
            # Allow a little branching: simply choose the branch with the greatest chain length
            new_visited = visited.copy()
            new_visited.add(nbr.GetIdx())
            branch_len = 1 + longest_aliphatic_chain_length(nbr, new_visited)
            if branch_len > max_len:
                max_len = branch_len
    return max_len

def is_fatty_amide(smiles: str):
    """
    Determines if a molecule is a fatty amide.
    
    We require that an amide bond (C(=O)N) is present such that
      (1) The acyl (carbonyl) side derives from a (possibly branched) aliphatic chain.
          We measure this by counting the maximum linear chain length starting from the acyl carbon (excluding the C=O oxygen)
          and require that length (from the carbon attached to the carbonyl) is at least acyl_min.
      (2) On the amide N, all substituents (other than the carbonyl connection) have a longest aliphatic chain length shorter than amine_max.
    
    We set acyl_min = 4 (so that even a branched isobutyl acyl chain as in dolichotheline qualifies)
    and amine_max = 4.
    
    Args:
      smiles (str): SMILES string of the molecule.
      
    Returns:
      bool: True if the molecule meets fatty amide criteria, False otherwise.
      str: Explanation of the decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # SMARTS for a basic amide bond (look for a carbonyl carbon sp2 bound to an oxygen by =, and bound to a nitrogen)
    amide_smarts = Chem.MolFromSmarts("C(=O)N")
    if amide_smarts is None:
        return False, "Error generating SMARTS pattern for amide"
        
    matches = mol.GetSubstructMatches(amide_smarts)
    if not matches:
        return False, "No amide bond (C(=O)N) found"
    
    # Define thresholds:
    acyl_min = 4    # minimum linear-chain length (from the chain attached to the carbonyl carbon)
    amine_max = 4   # maximum allowed chain length on substituents from the amide nitrogen
    
    # Process each amide match.
    # Each match is a tuple of indices: (carbonyl carbon, carbonyl oxygen, amide nitrogen)
    for match in matches:
        carbonyl = mol.GetAtomWithIdx(match[0])
        oxy = mol.GetAtomWithIdx(match[1])
        amide_nitrogen = mol.GetAtomWithIdx(match[2])
        
        # From the carbonyl carbon, determine the acyl side.
        # Exclude the carbonyl oxygen and the amide nitrogen.
        acyl_neighbors = [nbr for nbr in carbonyl.GetNeighbors()
                          if nbr.GetIdx() not in {oxy.GetIdx(), amide_nitrogen.GetIdx()}
                          and nbr.GetAtomicNum() == 6
                          and (not nbr.GetIsAromatic())
                          and (not nbr.IsInRing())]
        if not acyl_neighbors:
            continue  # this amide bond does not have an acyl substituent
        
        # For each acyl neighbor, compute the longest aliphatic chain length.
        acyl_chain_length = 0
        for nbr in acyl_neighbors:
            # starting count: include this neighbor
            chain_len = 1 + longest_aliphatic_chain_length(nbr, visited={carbonyl.GetIdx(), nbr.GetIdx()})
            if chain_len > acyl_chain_length:
                acyl_chain_length = chain_len
        
        # Check if the acyl chain meets the minimum threshold.
        if acyl_chain_length < acyl_min:
            # This amide does not have a sufficiently long acyl (fatty acid) group.
            continue
        
        # Now check the substituents on the amide nitrogen.
        # We want every substituent (aside from the carbonyl connection) to be modest (i.e. short aliphatic chain).
        valid_amine = True
        max_amine_chain = 0  # record the longest chain found on N substituents
        for nbr in amide_nitrogen.GetNeighbors():
            if nbr.GetIdx() == carbonyl.GetIdx():
                continue  # skip the carbonyl side
            if nbr.GetAtomicNum() == 6 and (not nbr.GetIsAromatic()) and (not nbr.IsInRing()):
                amine_chain = 1 + longest_aliphatic_chain_length(nbr, visited={amide_nitrogen.GetIdx(), nbr.GetIdx()})
                if amine_chain > max_amine_chain:
                    max_amine_chain = amine_chain
                if amine_chain >= amine_max:
                    valid_amine = False
                    break
        if not valid_amine:
            continue
        
        # If we reach here, we found a qualifying amide bond.
        return True, (f"Found fatty amide: acyl chain length {acyl_chain_length} "
                      f"and maximum N-side chain length {max_amine_chain}")
    
    return False, "Amide groups were found but none met fatty amide criteria (acyl chain too short or N-substituents too long)."


# Example usage (for testing):
if __name__ == "__main__":
    test_smiles = [
        # Should be classified as fatty amide:
        "CC(C)CC(=O)NCCc1c[nH]cn1",  # Dolichotheline 
        "CCCCCCCCC(=O)NCc1ccc(O)c(OC)c1",  # nonivamide
        "O=C(NCC(O)=O)CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CC",  # Docosahexaenoyl glycine
        # Others that might fail (or be rejected) by our criteria:
        "O=C(N[C@H](C(=O)N[C@@H](C(O)O)CCC(C)C)CCC(=O)N)[C@@H](NC(=O)C[C@H](O)CCCCCCCCCCC)CC(=O)N"  # a peptide-like molecule
    ]
    
    for smi in test_smiles:
        res, msg = is_fatty_amide(smi)
        print(f"SMILES: {smi}\nResult: {res}\nMessage: {msg}\n")