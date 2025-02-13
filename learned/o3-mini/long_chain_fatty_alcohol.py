"""
Classifies: CHEBI:17135 long-chain fatty alcohol
"""
#!/usr/bin/env python
"""
Classifies: long-chain fatty alcohol
Definition: A fatty alcohol with a chain length ranging from C13 to C22.
A valid long-chain fatty alcohol is expected to have an aliphatic –OH group attached to a non‐aromatic,
acyclic carbon that is part of the single longest (dominant) aliphatic chain. The longest chain must have
a length between 13 and 22 carbons and must account for most of the molecule's carbons.
Additionally, molecules that are carboxylic acids (or contain free –COOH groups) are disqualified.
"""

from rdkit import Chem

def is_long_chain_fatty_alcohol(smiles: str):
    """
    Determines if a molecule is a long-chain fatty alcohol.
    A long-chain fatty alcohol is defined as having an aliphatic –OH group attached to a carbon (which is in
    a dominant, contiguous aliphatic chain) and the longest such chain (composed solely of carbons that are
    not aromatic and not in any ring) has a length between 13 and 22. In addition, the dominant chain should
    represent the majority of the molecule's carbon atoms.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if the molecule is a long-chain fatty alcohol, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # --------- NEW CHECK: Disqualify molecules with a free carboxylic acid group ---------
    acid_query = Chem.MolFromSmarts("[CX3](=O)[OX1H]")
    if mol.HasSubstructMatch(acid_query):
        return False, "Molecule contains a free carboxylic acid group, so not a fatty alcohol"
    
    # Identify candidate alcohol groups (–OH). We use a SMARTS for [OX2H].
    alcohol_query = Chem.MolFromSmarts("[OX2H]")
    alcohol_matches = mol.GetSubstructMatches(alcohol_query)
    candidate_alcohol_carbons = set()
    if not alcohol_matches:
        return False, "No -OH functional group found"
    
    # For each -OH match, record the attached carbon (if non-aromatic)
    for match in alcohol_matches:
        o_idx = match[0]
        o_atom = mol.GetAtomWithIdx(o_idx)
        for nbr in o_atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6 and (not nbr.GetIsAromatic()):
                candidate_alcohol_carbons.add(nbr.GetIdx())
    
    if not candidate_alcohol_carbons:
        return False, "No aliphatic -OH group (attached to a non-aromatic carbon) found"
    
    # Determine the total number of carbon atoms
    total_carbons = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    if not total_carbons:
        return False, "No carbon atoms found"
    
    # Identify carbons that are in rings. (We want a contiguous chain in acyclic regions.)
    ring_info = mol.GetRingInfo()
    ring_carbons = set()
    for ring in ring_info.AtomRings():
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() == 6:
                ring_carbons.add(idx)
    
    # Build a set of carbon atoms that are both acyclic and non-aromatic.
    nonring_aliphatic = set()
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6 and (not atom.GetIsAromatic()) and atom.GetIdx() not in ring_carbons:
            nonring_aliphatic.add(atom.GetIdx())
    
    if not nonring_aliphatic:
        return False, "No acyclic aliphatic carbon chain found (all carbons are in rings or aromatic)"
    
    # Build a graph (dictionary) connecting the acyclic nonring aliphatic carbons.
    carbon_graph = {}
    for idx in nonring_aliphatic:
        atom = mol.GetAtomWithIdx(idx)
        neighbors = []
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6 and nbr.GetIdx() in nonring_aliphatic:
                neighbors.append(nbr.GetIdx())
        carbon_graph[idx] = neighbors

    # Find the longest simple path in the graph via DFS.
    longest_chain_length = 0
    longest_chain_path = []

    def dfs(node, visited):
        best_path = [node]
        for neighbor in carbon_graph.get(node, []):
            if neighbor not in visited:
                candidate_path = [node] + dfs(neighbor, visited | {neighbor})
                if len(candidate_path) > len(best_path):
                    best_path = candidate_path
        return best_path

    for idx in carbon_graph:
        path = dfs(idx, {idx})
        if len(path) > longest_chain_length:
            longest_chain_length = len(path)
            longest_chain_path = path

    if longest_chain_length == 0:
        return False, "Could not determine a contiguous aliphatic carbon chain"
    
    # Ensure that the longest chain length is within the desired range (C13 - C22)
    if not (13 <= longest_chain_length <= 22):
        return False, (f"Longest carbon chain length is {longest_chain_length}, "
                       f"which is not in the required range C13-C22")
    
    # Check if at least one candidate -OH (its attached carbon) is present in the longest chain.
    if not candidate_alcohol_carbons.intersection(longest_chain_path):
        return False, "The aliphatic -OH group is not attached to the dominant carbon chain"
    
    # Check that the dominant acyclic chain represents most of the molecule's carbons.
    total_carbon_count = len(total_carbons)
    ratio = longest_chain_length / total_carbon_count
    if ratio < 0.85:
        return False, (f"Longest chain (length {longest_chain_length}) accounts for only {ratio:.2f} "
                       f"of the total {total_carbon_count} carbons, so the molecule has significant additional fragments")
    
    return True, (f"Contains an aliphatic -OH attached to a dominant carbon chain of length "
                  f"{longest_chain_length} (which is {ratio:.2f} of all carbons) within C13-C22")

# Example usage:
# To test some SMILES strings, you can uncomment the following lines:
# print(is_long_chain_fatty_alcohol("CCCCCCCCCCCCCC(O)CCCC"))  # nonadecan-5-ol, expected True
# print(is_long_chain_fatty_alcohol("O1[C@@H]2[C@H](O)C(/C=C/[C@H](O)CCCCCCCC)=C([C@H]([C@H]12)O)CO"))  # Phomopoxide D, expected False (chain length 11)