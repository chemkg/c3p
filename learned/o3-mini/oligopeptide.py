"""
Classifies: CHEBI:25676 oligopeptide
"""
#!/usr/bin/env python
"""
Classifies: A peptide containing a relatively small number of amino acids (oligopeptide)

This improved approach tries to capture a contiguous peptide backbone.
Improvements over the previous method:
  1. We add an extra SMARTS pattern without stereochemical tags to capture peptide bonds
     that are not explicitly chiral.
  2. When “chaining” residues, we check that the bond linking a carbonyl in one residue
     to the amine in the next is not in a ring (which helps avoid false positives from
     beta‐lactam systems or similar cyclic amides).
We assume that an oligopeptide is a linear peptide with 2–10 amino acids.
"""

from rdkit import Chem

def is_oligopeptide(smiles: str):
    """
    Determines if a molecule is an oligopeptide based on its SMILES string.
    An oligopeptide is defined here as a peptide (linear, not cyclic)
    with at least 2 and no more than 10 amino acids.
    
    The algorithm is as follows:
      - Parse the SMILES.
      - Search for peptide backbone patterns. Instead of counting every "C(=O)N" (which
        may match nonpeptide rings), we look for substructures matching:
          •  N[C@@H](*)C(=O)
          •  N[C@H](*)C(=O)
          •  NC(*)C(=O)
      - We deduplicate matches using the index of the central (alpha) carbon.
      - We then build a connectivity (directed) graph among matches:
          An edge is drawn from residue A to residue B if the carbonyl carbon of A (position 2 in the match)
          is directly bonded to the amine nitrogen of B (position 0) AND that connecting bond is not in a ring.
      - We perform a DFS over this graph (avoiding cycles) to get the longest chain.
      - The number of residues is taken as (peptide bonds in the chain + 1). We require a chain
        length of 2–10 residues for an oligopeptide.
    
    Args:
      smiles (str): SMILES string of the molecule.
    
    Returns:
      bool: True if the molecule is an oligopeptide, False otherwise.
      str: Reason for the classification.
    """
    # Parse the molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define three SMARTS patterns; the third is less strict (no chiral annotation)
    patterns = [
        Chem.MolFromSmarts("N[C@@H](*)C(=O)"),
        Chem.MolFromSmarts("N[C@H](*)C(=O)"),
        Chem.MolFromSmarts("NC(*)C(=O)")
    ]
    
    # Collect all matches from all patterns.
    # Each match is a tuple: (N_index, alpha_C_index, carbonyl_C_index)
    all_matches = []
    for pat in patterns:
        if pat is None:
            continue
        ms = mol.GetSubstructMatches(pat)
        all_matches.extend(ms)
    
    if not all_matches:
        return False, "No peptide backbone patterns detected; not a peptide."
    
    # Deduplicate matches based on the central (alpha) carbon.
    # In a valid backbone match the alpha carbon should be unique.
    unique_matches = {}
    for match in all_matches:
        # match[1] is the index of the supposed alpha carbon.
        if match[1] not in unique_matches:
            unique_matches[match[1]] = match
    residue_matches = list(unique_matches.values())
    
    if not residue_matches:
        return False, "No unique peptide backbone residues detected."
    
    # Build a directed graph among residue matches.
    # For a pair of residues A and B, add an edge from A to B if:
    #   - The carbonyl carbon in A (match[2]) is directly bonded to the amine nitrogen in B (match[0])
    #   - AND that connecting bond is not in a ring.
    graph = {i: [] for i in range(len(residue_matches))}
    for i, match_i in enumerate(residue_matches):
        carbonyl_i = match_i[2]
        for j, match_j in enumerate(residue_matches):
            if i == j:
                continue
            amine_j = match_j[0]
            bond = mol.GetBondBetweenAtoms(carbonyl_i, amine_j)
            if bond and not bond.IsInRing():
                graph[i].append(j)
    
    # DFS to find the longest chain (each edge indicates one peptide bond).
    def dfs(node, visited):
        max_length = 0
        for neighbor in graph[node]:
            if neighbor in visited:
                continue
            length = 1 + dfs(neighbor, visited | {neighbor})
            if length > max_length:
                max_length = length
        return max_length

    longest_chain = 0
    for node in graph:
        chain_length = dfs(node, {node})
        if chain_length > longest_chain:
            longest_chain = chain_length

    # The number of residues is one more than the number of bonds in the chain.
    num_residues = longest_chain + 1

    if num_residues < 2:
        return False, "Only a single residue detected; not a peptide."
    if num_residues > 10:
        return False, f"Peptide has {num_residues} amino acids, which is too many to be considered an oligopeptide."
    
    return True, f"Peptide detected with {num_residues} amino acids; classified as an oligopeptide."

# Uncomment below to test an example, e.g., Arg-Arg-Phe:
# test_smiles = "O=C(N[C@@H](CC1=CC=CC=C1)C(O)=O)[C@@H](NC(=O)[C@@H](N)CCCN=C(N)N)CCCN=C(N)N"
# print(is_oligopeptide(test_smiles))