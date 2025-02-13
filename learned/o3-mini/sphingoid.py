"""
Classifies: CHEBI:35785 sphingoid
"""
#!/usr/bin/env python
"""
Classifies: Sphingoid
Definition: Sphinganine, its homologs and stereoisomers, and the hydroxy/unsaturated derivatives of these compounds.
The algorithm:
  1. Parse the SMILES and verify it is valid.
  2. Compute the longest contiguous chain of aliphatic (non‐ring) carbons. (A typical sphingoid has a long (≥14) alkyl chain.)
  3. Count the number of nitrogen atoms.
  4. Check whether any –OH group is present (if not, note that this may indicate a deoxy derivative).
  5. Look for a sphingoid “core” motif. Because sphingoid cores vary in exact linkage, several SMARTS patterns are tested:
       • A “full” motif: –OH–C–C(NH3/–NH2)–C–OH,
       • A deoxy variant missing one –OH,
       • An unsaturated variant,
       • And an “extended” version that allows an extra carbon spacer.
  6. If at least one of these patterns is found (and the chain length and nitrogen criteria met), return True.
  
If a molecule does not meet the requirements, a detailed reason is returned.
"""

from rdkit import Chem

# Helper: Compute longest contiguous chain of aliphatic (non‐ring, non‐aromatic) carbon atoms.
def _longest_aliphatic_chain_length(mol):
    eligible = set()
    for atom in mol.GetAtoms():
        # Only count carbons that are not in a ring and not aromatic.
        if atom.GetSymbol() == "C" and (not atom.IsInRing()) and (not atom.GetIsAromatic()):
            eligible.add(atom.GetIdx())
    # Build a neighbor dictionary for these eligible atoms.
    neighbors = {}
    for idx in eligible:
        atom = mol.GetAtomWithIdx(idx)
        neigh_idxs = []
        for neigh in atom.GetNeighbors():
            if neigh.GetIdx() in eligible:
                neigh_idxs.append(neigh.GetIdx())
        neighbors[idx] = neigh_idxs

    # DFS to determine longest connected path.
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
      - A long contiguous aliphatic chain (≥14 carbons among non‐ring carbons)
      - At least one nitrogen atom (as part of the head‐group)
      - The presence of one or more –OH groups in the molecule
      - The presence of a sphingoid core motif, as detected by several SMARTS patterns.
        Depending on the pattern a note is added whether it matches a full, deoxy, or unsaturated motif.
    
    Args:
      smiles (str): SMILES string of the molecule
    
    Returns:
      (bool, str): classification flag and a reason.
    """
    # 1. Parse the SMILES.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 2. Compute longest contiguous aliphatic (non‐ring) carbon chain.
    chain_len = _longest_aliphatic_chain_length(mol)
    if chain_len < 14:
        return False, f"Longest continuous aliphatic carbon chain is only {chain_len} (need at least 14)"

    # 3. Count nitrogen atoms.
    nitrogen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    if nitrogen_count < 1:
        return False, "No nitrogen found (expected an amino or acylamino group in sphingoid compounds)"

    # 4. Check for hydroxyl groups.
    oh_pattern = Chem.MolFromSmarts("[OX2H]")
    has_oh = mol.HasSubstructMatch(oh_pattern)
    deoxy_note = "" if has_oh else " (no OH found, so may be a deoxy derivative)"

    # 5. Look for sphingoid core motif.
    # We define several SMARTS patterns to capture the variation in sphingoid headgroups.
    patterns = [
       ("full", "[OX2H][C;!R]-[C;!R]([NX3])[C;!R][OX2H]"),
       ("deoxy", "[C;!R]-[C;!R]([NX3])[C;!R][OX2H]"),
       ("unsaturated", "[C;!R]=[C;!R]([NX3])"),
       ("extended", "[OX2H][C;!R]-[C;!R]-[C;!R]([NX3])[CH2][OX2H]")
    ]
    match_type = None
    for name, smarts in patterns:
        patt = Chem.MolFromSmarts(smarts)
        if patt is None:
            continue  # skip if pattern did not compile
        if mol.HasSubstructMatch(patt):
            match_type = name
            break

    if not match_type:
        return False, "No sphingoid core motif found in molecule"

    # Interpret the pattern match.
    if match_type == "full" or match_type == "extended":
        core_desc = "Matches full sphingoid headgroup motif"
    elif match_type == "deoxy":
        core_desc = "Matches deoxy sphingoid motif (lacking one OH in headgroup)"
    elif match_type == "unsaturated":
        core_desc = "Matches unsaturated sphingoid motif"
    else:
        core_desc = "Matches sphingoid core motif"

    reason = (f"{core_desc}{deoxy_note}. Long alkyl chain of {chain_len} carbons "
              f"and {nitrogen_count} nitrogen(s) detected.")
    return True, reason

# Example usage and testing with several SMILES from the list:
if __name__ == "__main__":
    examples = [
        ("CCCCCCCCCCCCCCCCCC(=O)N[C@@H](CO[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O)[C@H](O)[C@H](O)CCCCCCCCCCC(C)C",
         "N-octadecanoyl-1-O-beta-D-glucosyl-4-hydroxy-15-methylhexadecasphinganine"),
        ("CCCCCCCCCCCC[C@@H](O)[C@@H](O)[C@@H]([NH3+])CO",
         "C16 phytosphingosine(1+)"),
        ("CCCCCCCCCCCC\\C=C\\[C@@H](O)[C@@H](N)CO",
         "tetradecasphingosine"),
        ("CCCCCCCCCCCCCCCC(=O)CN", "1-deoxymethyl-3-dehydrosphinganine"),
        ("CC(C)CCCCCCCCCCCC(=O)[C@@H](N)CO", "3-dehydro-15-methylhexadecasphinganine")
    ]
    
    for smi, name in examples:
        flag, reason = is_sphingoid(smi)
        print(f"Name: {name}\nSMILES: {smi}\n  --> Sphingoid: {flag}, Reason: {reason}\n")