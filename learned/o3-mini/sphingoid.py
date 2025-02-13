"""
Classifies: CHEBI:35785 sphingoid
"""
#!/usr/bin/env python
"""
Classifies: Sphingoid
Definition: Sphinganine, its homologs and stereoisomers, 
         and the hydroxy and unsaturated derivatives of these compounds.
         
The algorithm:
  1. Parse the SMILES; if invalid, return an error.
  2. Compute the longest contiguous chain of aliphatic (non‐ring) carbons.
     (A true sphingoid (sphinganine, sphingosine, phytosphingosine etc.)
      normally has a long (>~14) alkyl chain.)
  3. Check that there is at least one nitrogen.
  4. Check for the presence of an –OH group.
  5. Look for a polar “head‐group” motif by scanning for one of several 
     SMARTS patterns that roughly correspond to the sphingoid core (e.g.
     HO–CH–CH(NH2)–CH2OH or derivatives).
  6. If any of the core motifs are found and the chain plus N requirements 
     are met, classify as a sphingoid.
     
Since sphingoid compounds may be “deoxy” (lack one OH) or unsaturated, the code 
allows different substructure patterns and prints a note if no –OH is found.
If the criteria are not all met, the function returns False with the reason.
"""

from rdkit import Chem

# Helper: find the longest contiguous chain of aliphatic (non‐ring) carbons.
# We define eligible carbons as those with symbol "C" that are not aromatic and not in a ring.
def _longest_aliphatic_chain_length(mol):
    # Gather indices of eligible carbon atoms.
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
    
    The function requires that:
      - The molecule has a long contiguous aliphatic chain (≥14 carbons out of non‐ring carbons);
      - It contains at least one nitrogen (to allow an amino or acylamino group);
      - It has at least one hydroxyl –OH (or else it is noted as possibly deoxy);
      - It contains a “core” polar fragment (one of several SMARTS patterns) that roughly matches
        the sphingoid head group (for example, a motif like HO–CH–CH(NH2)–CH2OH, or variants).
    
    Args:
      smiles (str): a SMILES string
    
    Returns:
      (bool, str): Classification flag plus a reason.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Step 1: Verify long aliphatic chain.
    chain_len = _longest_aliphatic_chain_length(mol)
    if chain_len < 14:
        return False, f"Longest continuous aliphatic carbon chain is only {chain_len} (need at least 14)"
    
    # Step 2: Check for nitrogen presence.
    nitro_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    if nitro_count < 1:
        return False, "No nitrogen found (expected an amino or acylamino group in sphingoid compounds)"
    
    # Step 3: Check for -OH groups.
    oh_pattern = Chem.MolFromSmarts("[OX2H]")
    has_oh = mol.HasSubstructMatch(oh_pattern)
    deoxy_note = "" if has_oh else " (no OH found, so may be a deoxy derivative)"
    
    # Step 4: Look for a sphingoid polar "head‐group" motif.
    # We define three alternative SMARTS:
    # Pattern A: Full headgroup (expects a hydroxyl at both ends of a 3–4 carbon unit)
    pattern_A = Chem.MolFromSmarts("[OX2H][C;!R]-[C;!R]([NX3;H2,H1])- [C;!R][OX2H]")
    # Pattern B: A deoxy variant lacking the first OH.
    pattern_B = Chem.MolFromSmarts("[C;!R]-[C;!R]([NX3;H2,H1])-[C;!R][OX2H]")
    # Pattern C: An unsaturated variant where a C=C is present before the amino carbon.
    pattern_C = Chem.MolFromSmarts("[C;!R]=[C;!R]([NX3;H2,H1])")
    # (Remove extraneous spaces in SMARTS)
    # Note: SMARTS can be finicky; these are heuristic patterns and might miss borderline cases.
    
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
    
    # If all conditions are met, classify as sphingoid.
    reason = f"{core_match}{deoxy_note}. Long alkyl chain of {chain_len} carbons and {nitro_count} nitrogen(s) detected."
    return True, reason

# Example usage:
if __name__ == "__main__":
    examples = [
        # True positives (from provided examples)
        ("CCCCCCCCCCCCCCCCCC(=O)N[C@@H](CO[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O)[C@H](O)[C@H](O)CCCCCCCCCCC(C)C",
         "N-octadecanoyl-1-O-beta-D-glucosyl-4-hydroxy-15-methylhexadecasphinganine"),
        ("CCCCCCCCCCCC[C@@H](O)[C@@H](O)[C@@H]([NH3+])CO",
         "C16 phytosphingosine(1+)"),
        ("CCCCCCCCCCCC\\C=C\\[C@@H](O)[C@@H](N)CO",
         "tetradecasphingosine"),
        # A false negative candidate from report:
        ("CCCCCCCCCCCCCCCC(=O)CN", "1-deoxymethyl-3-dehydrosphinganine"),
    ]
    
    for smi, name in examples:
        flag, reason = is_sphingoid(smi)
        print(f"Name: {name}\nSMILES: {smi}\n --> Sphingoid: {flag}, Reason: {reason}\n")