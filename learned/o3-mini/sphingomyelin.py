"""
Classifies: CHEBI:64583 sphingomyelin
"""
"""
Classifies: Sphingomyelin

Definition:
  Sphingomyelins are a class of phospholipids in which the amino group of a 
  sphingoid base is in amide linkage with a fatty acid, while the terminal 
  hydroxyl group of the sphingoid base is esterified to phosphorylcholine.
  
This implementation first looks for a core motif – a fatty acid carbonyl ([CX3](=O))
directly linked (via an amide bond) to a chiral carbon which is substituted with an 
attachment to a phosphorylcholine head group. In order to reduce false positives, 
the code then follows the bond connectivity to check that the “fatty acyl” branch 
(from the carbonyl carbon) and the “sphingoid chain” (from the chiral carbon) each 
have a sufficiently long, contiguous chain of carbon atoms. 

We allow for slight variation in the phosphate representation (neutral vs charged)
and for the chiral annotation ([C@H] or [C@@H]).
"""

from rdkit import Chem

# Helper function to compute the length of a carbon chain (number of carbon atoms 
# in the longest continuous path) that starts at a given atom. We follow bonds between
# carbons only and “avoid” going back to the atom we came from.
def chain_length(atom, coming_from_idx, mol, visited=None):
    if visited is None:
        visited = set()
    visited.add(atom.GetIdx())
    max_len = 1  # count the current atom
    for nbr in atom.GetNeighbors():
        if nbr.GetIdx() == coming_from_idx:
            continue
        # Only follow carbons
        if nbr.GetAtomicNum() != 6:
            continue
        if nbr.GetIdx() in visited:
            continue
        # We follow regardless of bond order (double bonds in unsaturated chains are allowed)
        current_len = 1 + chain_length(nbr, atom.GetIdx(), mol, visited.copy())
        if current_len > max_len:
            max_len = current_len
    return max_len

def is_sphingomyelin(smiles: str):
    """
    Determines if a molecule is a sphingomyelin based on its SMILES string.
    
    The strategy is as follows:
      1. Parse the SMILES string with RDKit.
      2. Look for the core sphingomyelin motif using SMARTS.
         The SMARTS requires a fatty acid carbonyl ([CX3](=O)) attached via an amide (N)
         to a chiral carbon ([C@H] or [C@@H]) that carries a phosphorylcholine head group.
         We allow two variants for the phosphate group (neutral or charged).
      3. For any match, we examine two “chain” lengths:
           a) From the carbonyl carbon, we follow the fatty acyl chain (excluding the amide N),
              and require a minimum chain length (here set to ≥8 carbons).
           b) From the chiral sphingoid carbon (the one directly attached to the amide and 
              that bears the phosphate substituent) we follow the “sphingoid” chain – that is,
              we discard the branch containing the phosphorylcholine (attached via an oxygen)
              and select the remaining carbon branch. This chain must also be at least 8 carbons long.
      4. If both chain lengths meet the threshold and the core motif was found, the molecule
         is classified as sphingomyelin.
    
    Args:
      smiles (str): SMILES string of the molecule.
      
    Returns:
      bool: True if the molecule is classified as sphingomyelin, False otherwise.
      str: Reason for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define variants of the sphingomyelin core SMARTS:
    core_patterns = [
        "[CX3](=O)N[C@H](COP(=O)(O)OCC[N+](C)(C)C)",
        "[CX3](=O)N[C@@H](COP(=O)(O)OCC[N+](C)(C)C)",
        "[CX3](=O)N[C@H](COP([O-])(=O)OCC[N+](C)(C)C)",
        "[CX3](=O)N[C@@H](COP([O-])(=O)OCC[N+](C)(C)C)"
    ]
    
    # Minimum required chain length (number of carbons) for each branch.
    MIN_CHAIN_LENGTH = 8
    
    valid_core_found = False
    reasons = []
    for smarts in core_patterns:
        core_query = Chem.MolFromSmarts(smarts)
        if core_query is None:
            continue
        matches = mol.GetSubstructMatches(core_query)
        if not matches:
            continue
        
        # We iterate through each match variant.
        for match in matches:
            # For our SMARTS the assumed mapping is:
            #   match[0] = carbonyl carbon (fatty acid C=O)
            #   match[1] = amide nitrogen
            #   match[2] = chiral carbon on sphingoid base (carries head group)
            #   match[3] = oxygen that attaches the phosphate group
            if len(match) < 4:
                continue
            carbonyl_idx = match[0]
            amideN_idx   = match[1]
            chiral_idx   = match[2]
            oxy_idx      = match[3]  # oxygen linking to phosphate head
            
            atom_carbonyl = mol.GetAtomWithIdx(carbonyl_idx)
            atom_chiral   = mol.GetAtomWithIdx(chiral_idx)
            
            # --- 1. Check fatty acyl chain length:
            # From the carbonyl carbon, select the neighbor that is not the amide nitrogen.
            fatty_chain_lengths = []
            for nbr in atom_carbonyl.GetNeighbors():
                if nbr.GetIdx() == amideN_idx:
                    continue
                if nbr.GetAtomicNum() == 6:
                    length = chain_length(nbr, carbonyl_idx, mol)
                    fatty_chain_lengths.append(length)
            if not fatty_chain_lengths:
                reasons.append("Core motif found but unable to identify fatty acyl chain from carbonyl carbon")
                continue
            fatty_chain_max = max(fatty_chain_lengths)
            
            # --- 2. Check sphingoid chain length:
            # From the chiral carbon (which bears the phosphorylcholine via oxy_idx),
            # choose neighbors that are carbons and are not the amide N and not the phosphate oxygen.
            sphingo_chain_lengths = []
            for nbr in atom_chiral.GetNeighbors():
                if nbr.GetIdx() in (amideN_idx, oxy_idx):
                    continue
                if nbr.GetAtomicNum() == 6:
                    length = chain_length(nbr, chiral_idx, mol)
                    sphingo_chain_lengths.append(length)
            if not sphingo_chain_lengths:
                reasons.append("Core motif found but unable to identify sphingoid chain from chiral carbon")
                continue
            sphingo_chain_max = max(sphingo_chain_lengths)
            
            # Check if both chain lengths are above thresholds.
            if fatty_chain_max < MIN_CHAIN_LENGTH:
                reasons.append(f"Fatty acyl chain is too short (found {fatty_chain_max} carbons, need ≥{MIN_CHAIN_LENGTH})")
                continue
            if sphingo_chain_max < MIN_CHAIN_LENGTH:
                reasons.append(f"Sphingoid chain is too short (found {sphingo_chain_max} carbons, need ≥{MIN_CHAIN_LENGTH})")
                continue
            
            # If we reach here, we have a core match with valid chain lengths.
            valid_core_found = True
            return True, ("Contains sphingomyelin core with valid fatty acyl and sphingoid chains "
                          f"(fatty chain length: {fatty_chain_max} C; sphingoid chain length: {sphingo_chain_max} C)")
    
    if not valid_core_found:
        # If some matches were found but failed chain‐length tests, list one of the reasons.
        if reasons:
            return False, "Core motif found but " + "; ".join(reasons)
        else:
            return False, ("Missing sphingomyelin core motif: expected a fatty acid carbonyl (C(=O)) directly bound via an "
                           "amide linkage to a chiral carbon carrying a phosphorylcholine head group")
    
    # Fallback (should not normally occur)
    return False, "Unable to classify molecule as sphingomyelin"

# Example usage (uncomment for testing):
# test_smiles = "CCCCCCCCCCCCCCCCCCCC(=O)N[C@@H](COP([O-])(=O)OCC[N+](C)(C)C)[C@H](O)[C@H](O)CCCCCCCCCCC(C)C"
# result, reason = is_sphingomyelin(test_smiles)
# print(result, reason)