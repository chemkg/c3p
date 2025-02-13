"""
Classifies: CHEBI:25676 oligopeptide
"""
"""
Classifies: A peptide containing a relatively small number of amino acids (oligopeptide)

In this improved approach we attempt to detect a contiguous peptide backbone
by matching a more specific SMARTS pattern (one that requires an α‐chiral center)
for internal residues. Then, by “chaining” together matches that share a common
amide bond we compute a peptide chain length. We assume a peptide with 2–10 residues
is an oligopeptide.
"""

from rdkit import Chem

def is_oligopeptide(smiles: str):
    """
    Determines if a molecule is an oligopeptide based on its SMILES string.
    An oligopeptide is defined as a (linear) peptide with at least 2 but no more than 10 amino acids.
    
    The determination here is based on identifying a peptide backbone pattern.
    Instead of counting every "C(=O)N" (which may give false positives), we look for
    residues that show an amine attached to a chiral α‐carbon, which in turn bears a carbonyl.
    Specifically, we search for substructures matching either "N[C@@H](*)C(=O)" or
    "N[C@H](*)C(=O)". These patterns are expected to occur in a peptide backbone.
    
    We then build a connectivity graph among the residue matches,
    “linking” two residues if the carbonyl carbon of one is directly bonded to the amine
    nitrogen of the next. The longest resulting chain’s length, plus one for the N-terminal residue,
    is taken as the number of amino acids. If that number is between 2 and 10 (inclusive) we classify
    the molecule as an oligopeptide.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule is an oligopeptide, False otherwise.
        str: Explanation for the classification.
    """
    # Parse the SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define two SMARTS patterns that indicate an internal residue in a peptide backbone.
    # They require an amine bound to a chiral (alpha) carbon that is attached to a C=O.
    pattern1 = Chem.MolFromSmarts("N[C@@H](*)C(=O)")
    pattern2 = Chem.MolFromSmarts("N[C@H](*)C(=O)")
    
    # Get all matches from both patterns
    matches1 = mol.GetSubstructMatches(pattern1)
    matches2 = mol.GetSubstructMatches(pattern2)
    # Combine and deduplicate (each match is a tuple of atom indices: (N, C_alpha, C_carbonyl))
    all_matches = list({m for m in matches1 + matches2})
    if not all_matches:
        return False, "No peptide backbone patterns detected; not a peptide."
    
    # Build a directed graph among the residue matches.
    # We will consider an edge from match A to match B if the carbonyl carbon (atom index #2) of A
    # is directly bonded (via a single bond) to the amine nitrogen (atom index #0) of B.
    graph = {i: [] for i in range(len(all_matches))}
    for i, match_i in enumerate(all_matches):
        carbonyl_i = match_i[2]
        for j, match_j in enumerate(all_matches):
            if i == j:
                continue
            amine_j = match_j[0]
            # Check if there is a bond between carbonyl_i and amine_j:
            if mol.GetBondBetweenAtoms(carbonyl_i, amine_j):
                graph[i].append(j)
    
    # Now, find the length of the longest path through this graph (each edge indicates one peptide bond).
    # A simple DFS (depth-first search) is used. We do not allow cycles.
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

    # The number of peptide bonds detected in the chain is 'longest_chain'
    # The number of amino acid residues is bonds + 1.
    num_residues = longest_chain + 1

    # We require at least 2 residues for a peptide (an amino acid alone is not one)
    if num_residues < 2:
        return False, "Only a single residue detected; not a peptide."
    # We set the threshold for an oligopeptide as no more than 10 amino acids.
    if num_residues > 10:
        return False, f"Peptide has {num_residues} amino acids, which is too many to be considered an oligopeptide."
    
    return True, f"Peptide detected with {num_residues} amino acids; classified as an oligopeptide."


# Uncomment below to run a test example (e.g. Arg-Arg-Phe)
# test_smiles = "O=C(N[C@@H](CC1=CC=CC=C1)C(O)=O)[C@@H](NC(=O)[C@@H](N)CCCN=C(N)N)CCCN=C(N)N"
# print(is_oligopeptide(test_smiles))