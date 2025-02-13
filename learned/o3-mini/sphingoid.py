"""
Classifies: CHEBI:35785 sphingoid
"""
#!/usr/bin/env python
"""
Classifies: Sphingoid
Definition: Sphinganine, its homologs and stereoisomers, and the hydroxy and unsaturated derivatives of these compounds.
The algorithm:
  1. Parse the SMILES; if invalid, return an error.
  2. Compute the longest contiguous chain of aliphatic (non‐ring) carbons.
     (A true sphingoid normally has a long (≥14) alkyl chain.)
  3. Check that there is at least one nitrogen.
  4. Check for the presence of an –OH group.
  5. Look for a polar “head‐group” motif by scanning for one of several SMARTS patterns that roughly correspond 
     to the sphingoid core (for example: HO–CH–CH(NH2)–CH2OH or slightly modified variants).
  6. If any of the core motifs are found and the chain plus nitrogen requirements met, classify as sphingoid.
     
If the –OH is absent, a note is added that the compound may be a deoxy derivative.
"""

from rdkit import Chem

# Helper: Find the longest contiguous chain of aliphatic (non‐ring) carbons.
def _longest_aliphatic_chain_length(mol):
    # Gather indices of eligible carbon atoms (non-ring, non-aromatic).
    eligible = set()
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == "C" and (not atom.IsInRing()) and (not atom.GetIsAromatic()):
            eligible.add(atom.GetIdx())
    
    # Build a neighbor dictionary for eligible carbons.
    neighbors = {}
    for idx in eligible:
        atom = mol.GetAtomWithIdx(idx)
        neigh_idxs = []
        for neigh in atom.GetNeighbors():
            if neigh.GetIdx() in eligible:
                neigh_idxs.append(neigh.GetIdx())
        neighbors[idx] = neigh_idxs

    # Use DFS to compute the longest path in the undirected graph.
    def dfs(current, visited):
        max_len = 1
        for nb in neighbors.get(current, []):
            if nb not in visited:
                path_len = 1 + dfs(nb, visited | {nb})
                if path_len > max_len:
                    max_len = path_len
        return max_len

    longest = 0
    for idx in eligible:
        chain_len = dfs(idx, {idx})
        if chain_len > longest:
            longest = chain_len
    return longest

def is_sphingoid(smiles: str):
    """
    Determines whether a molecule belongs to the sphingoid class 
    (sphinganine, its homologs and stereoisomers, and the hydroxy/unsaturated derivatives).
    
    Requirements:
      - Long contiguous aliphatic chain (≥14 carbons among non‐ring carbons)
      - At least one nitrogen (for the amino or acylamino group)
      - Presence (or note the absence) of an –OH group 
      - At least one "core" polar motif corresponding to the sphingoid headgroup.
    
    Args:
      smiles (str): SMILES string of the molecule
    
    Returns:
      (bool, str): classification flag and a reason.
    """
    # Parse the SMILES.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Step 1: Compute the longest contiguous aliphatic chain length.
    chain_len = _longest_aliphatic_chain_length(mol)
    if chain_len < 14:
        return False, f"Longest continuous aliphatic carbon chain is only {chain_len} (need at least 14)"
    
    # Step 2: Count the nitrogen atoms.
    nitrogen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    if nitrogen_count < 1:
        return False, "No nitrogen found (expected an amino or acylamino group in sphingoid compounds)"
    
    # Step 3: Check for hydroxyl groups.
    oh_pattern = Chem.MolFromSmarts("[OX2H]")
    has_oh = mol.HasSubstructMatch(oh_pattern)
    deoxy_note = "" if has_oh else " (no OH found, so may be a deoxy derivative)"
    
    # Step 4: Look for a sphingoid polar headgroup motif.
    # Patterns correspond roughly to HO–CH–CH(NH2)–CH2OH or slight variations.
    # Corrected SMARTS (removed extraneous spaces):
    pattern_A = Chem.MolFromSmarts("[OX2H][C;!R]-[C;!R]([NX3;H2,H1])[C;!R][OX2H]")
    pattern_B = Chem.MolFromSmarts("[C;!R]-[C;!R]([NX3;H2,H1])[C;!R][OX2H]")
    pattern_C = Chem.MolFromSmarts("[C;!R]=[C;!R]([NX3;H2,H1])")
    
    # Check that SMARTS compiled properly.
    if pattern_A is None or pattern_B is None or pattern_C is None:
        return False, "SMARTS pattern compilation failed"
    
    matches_A = mol.GetSubstructMatches(pattern_A)
    matches_B = mol.GetSubstructMatches(pattern_B)
    matches_C = mol.GetSubstructMatches(pattern_C)
    
    core_match = None
    if matches_A:
        core_match = "Matches full sphingoid headgroup motif"
    elif matches_B:
        core_match = "Matches deoxy sphingoid motif (lacking one OH in headgroup)"
    elif matches_C:
        core_match = "Matches unsaturated sphingoid motif"
    
    if not core_match:
        return False, "No sphingoid core motif found in molecule"
    
    reason = (f"{core_match}{deoxy_note}. Long alkyl chain of {chain_len} carbons "
              f"and {nitrogen_count} nitrogen(s) detected.")
    return True, reason

# Example usage:
if __name__ == "__main__":
    # List a few example SMILES strings from the provided list.
    examples = [
        ("CCCCCCCCCCCCCCCCCC(=O)N[C@@H](CO[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O)[C@H](O)[C@H](O)CCCCCCCCCCC(C)C",
         "N-octadecanoyl-1-O-beta-D-glucosyl-4-hydroxy-15-methylhexadecasphinganine"),
        ("CCCCCCCCCCCC[C@@H](O)[C@@H](O)[C@@H]([NH3+])CO",
         "C16 phytosphingosine(1+)"),
        ("CCCCCCCCCCCC\\C=C\\[C@@H](O)[C@@H](N)CO",
         "tetradecasphingosine"),
        ("CCCCCCCCCCCCCCCC(=O)CN", "1-deoxymethyl-3-dehydrosphinganine") # a candidate that may fail headgroup match
    ]
    
    for smi, name in examples:
        flag, reason = is_sphingoid(smi)
        print(f"Name: {name}\nSMILES: {smi}\n  --> Sphingoid: {flag}, Reason: {reason}\n")