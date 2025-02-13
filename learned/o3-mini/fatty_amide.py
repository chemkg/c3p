"""
Classifies: CHEBI:29348 fatty amide
"""
"""
Classifies: Fatty Amide
Definition:
  A monocarboxylic acid amide derived from a fatty acid.
  Here we require that the amide bond C(=O)N has an acyl (carbonyl) side
  having at least a minimum linear chain length (acyl_min) of 4 (allowing branched chains)
  and that each substituent on the amide nitrogen (other than the carbonyl link)
  is comparatively short – with its longest aliphatic chain being no longer than amine_max (4).
  
Note:
  This heuristic was developed after earlier attempts showed many false positives (eg, from peptides)
  and false negatives (eg, missing dolichotheline). In our revised approach we also allow an N-chain of
  length exactly 4, but reject if any is longer.
"""

from rdkit import Chem

def longest_aliphatic_chain_length(atom, visited):
    """
    Recursively computes the longest linear chain length starting from the given atom.
    We only traverse carbon atoms that are:
      - nonaromatic,
      - not in any ring.
    visited is a set of atom indices already visited.
    The count always includes the starting atom.
    """
    max_len = 1  # count the current atom
    for nbr in atom.GetNeighbors():
        if nbr.GetAtomicNum() == 6 and (not nbr.GetIsAromatic()) and (not nbr.IsInRing()):
            if nbr.GetIdx() in visited:
                continue
            new_visited = visited.copy()
            new_visited.add(nbr.GetIdx())
            branch_len = 1 + longest_aliphatic_chain_length(nbr, new_visited)
            if branch_len > max_len:
                max_len = branch_len
    return max_len

def is_fatty_amide(smiles: str):
    """
    Determines if a molecule is a fatty amide.
    
    For each amide bond (C(=O)N) detected, the following is required:
      (1) On the carbonyl (acyl) side: at least one substituent that is an aliphatic chain
          having a linear-chain length (starting from the atom bound to the carbonyl carbon)
          of at least acyl_min.
      (2) On the amide nitrogen: every substituent (other than the carbonyl connection)
          has a longest aliphatic chain length of at most amine_max.
          
    Thresholds:
      acyl_min = 4   (so even branched isobutyl chains qualify)
      amine_max = 4  (N–substituents must be short; a chain length of 4 is tolerated, >4 is not)
    
    Args:
      smiles (str): SMILES string of the molecule.
      
    Returns:
      bool: True if at least one amide bond meets fatty amide criteria, False otherwise.
      str: Explanation of the decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Pattern for a basic amide bond: carbonyl carbon bound to an oxygen with a double bond and to a nitrogen.
    amide_smarts = Chem.MolFromSmarts("C(=O)N")
    if amide_smarts is None:
        return False, "Error generating SMARTS pattern for amide"
        
    matches = mol.GetSubstructMatches(amide_smarts)
    if not matches:
        return False, "No amide bond (C(=O)N) found"
    
    # Define thresholds:
    acyl_min = 4    # minimum chain length from the carbon directly bonded to the carbonyl
    amine_max = 4   # maximum allowed chain length on each substituent from the amide nitrogen
    
    # Loop over each amide bond match; a match returns (carbonylC, carbonylO, amideN)
    for match in matches:
        carbonyl = mol.GetAtomWithIdx(match[0])
        oxy = mol.GetAtomWithIdx(match[1])
        amide_nitrogen = mol.GetAtomWithIdx(match[2])
        
        # --- Process the acyl side ---
        # From the carbonyl, get the substituents other than the oxygen (of the C=O) and the amide N.
        acyl_neighbors = [nbr for nbr in carbonyl.GetNeighbors()
                          if nbr.GetIdx() not in {oxy.GetIdx(), amide_nitrogen.GetIdx()}
                          and nbr.GetAtomicNum() == 6
                          and (not nbr.GetIsAromatic())
                          and (not nbr.IsInRing())]
        if not acyl_neighbors:
            # No proper acyl substituent found for this amide bond.
            continue
        
        acyl_chain_length = 0
        for nbr in acyl_neighbors:
            # Count the chain beginning at this neighbor.
            chain_len = 1 + longest_aliphatic_chain_length(nbr, visited={carbonyl.GetIdx(), nbr.GetIdx()})
            if chain_len > acyl_chain_length:
                acyl_chain_length = chain_len
        
        if acyl_chain_length < acyl_min:
            # This amide bond does not have a long enough fatty acyl chain.
            continue
        
        # --- Process the substituents on the amide nitrogen ---
        valid_amine = True
        max_amine_chain = 0  # record the longest chain on N substituents
        for nbr in amide_nitrogen.GetNeighbors():
            if nbr.GetIdx() == carbonyl.GetIdx():
                continue  # skip the connection back to the carbonyl
            if nbr.GetAtomicNum() == 6 and (not nbr.GetIsAromatic()) and (not nbr.IsInRing()):
                # Get the chain length starting from this substituent.
                chain_len = 1 + longest_aliphatic_chain_length(nbr, visited={amide_nitrogen.GetIdx(), nbr.GetIdx()})
                if chain_len > max_amine_chain:
                    max_amine_chain = chain_len
                # Allow chain length equal to amine_max; reject only if it is longer.
                if chain_len > amine_max:
                    valid_amine = False
                    break
        if not valid_amine:
            continue
        
        # If found one amide bond that meets criteria, return True.
        return True, (f"Found fatty amide: acyl chain length {acyl_chain_length} "
                      f"and maximum N-side chain length {max_amine_chain}")
    
    return False, "Amide groups were found but none met fatty amide criteria (acyl chain too short or N-substituents too long)."


# Example usage (for testing):
if __name__ == "__main__":
    test_smiles = [
        # True positives examples:
        "CC(C)CC(=O)NCCc1c[nH]cn1",  # Dolichotheline: acyl chain length 4; N-side ~3 (accepted)
        "CCCCCCCCC(=O)NCc1ccc(O)c(OC)c1",  # nonivamide: acyl chain length 9; N-side ~2
        "O=C(NCC(O)=O)CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CC",  # Docosahexaenoyl glycine
        "C(CCC)C[C@@H](/C=C/C=C\\C/C=C\\C/C=C\\CCCC(N[C@H](C(O)=O)C)=O)OO",  # N-[(15S)-hydroperoxy-...]-alanine
        # Example of a false positive from previous attempts (a peptide) should now be rejected:
        "OC(C(N)C(=O)NC(CC1=CC=C(O)C=C1)C(=O)NCC(O)=O)C",  # A dipeptide-like structure
        # Example false negative from previous attempt:
        "C(=C/[C@H](C/C=C\\CCCCC)OO)\\C=C/C/C=C\\CCCC(=O)NCCCC(O)=O",  
          # N-[(12S)-hydroperoxy-(5Z,8Z,10E,14Z)-icosatetraenoyl]-gamma-aminobutanoic acid
    ]
    
    for smi in test_smiles:
        result, reason = is_fatty_amide(smi)
        print(f"SMILES: {smi}\nResult: {result}\nReason: {reason}\n")