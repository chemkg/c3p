"""
Classifies: CHEBI:59644 oxo fatty acid
"""
"""
Classifies: Oxo Fatty Acid
Definition: Any fatty acid containing at least one aldehydic or ketonic (oxo) group, in addition to the carboxylic acid group.
The molecule should have a terminal (free) carboxylic acid group and exhibit a long, predominantly acyclic alkyl chain.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_oxo_fatty_acid(smiles: str):
    """
    Determines if a molecule is an oxo fatty acid based on its SMILES string.
    An oxo fatty acid, in this context, is defined as a long-chain (predominantly acyclic) fatty acid that contains
    at least one additional aldehyde or ketone group outside of a terminal carboxylic acid functionality.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule meets the oxo fatty acid criteria, False otherwise.
        str: A reason explaining the classification decision.
    """

    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # ---------------------------------------------
    # (A) Check for a carboxylic acid group.
    # We require a free carboxylic acid group. The standard SMARTS pattern for carboxylic acid is:
    # [CX3](=O)[OX2H1]
    acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    acid_matches = mol.GetSubstructMatches(acid_pattern)
    if not acid_matches:
        return False, "No carboxylic acid group detected; not a fatty acid"

    # For each acid group match, check if the acid carbon is terminal.
    # We get the set of indices that are carboxyl carbons which are terminal, i.e. they have only one carbon neighbor.
    terminal_acid_indices = set()
    for match in acid_matches:
        acid_carbon = match[0]  # first atom in the SMARTS: the carbon of the acid
        # Count neighboring atoms which are carbon (atomic number 6);
        # ignore oxygen neighbors.
        carbon_neighbor_count = 0
        for nbr in mol.GetAtomWithIdx(acid_carbon).GetNeighbors():
            if nbr.GetAtomicNum() == 6:
                carbon_neighbor_count += 1
        # A terminal carboxylic acid should have exactly one carbon neighbor.
        if carbon_neighbor_count == 1:
            terminal_acid_indices.add(acid_carbon)
            
    if not terminal_acid_indices:
        return False, "Carboxylic acid group is not terminal; fatty acids require a free (terminal) acid group"
    
    # ---------------------------------------------
    # (B) Identify additional oxo groups (aldehyde or ketone) that lie outside of the acid.
    # Define SMARTS patterns for ketone and aldehyde.
    ketone_pattern = Chem.MolFromSmarts("[#6][CX3](=O)[#6]")
    aldehyde_pattern = Chem.MolFromSmarts("[#6][CX3H](=O)")
    
    ketone_matches = mol.GetSubstructMatches(ketone_pattern)
    aldehyde_matches = mol.GetSubstructMatches(aldehyde_pattern)
    
    # Check that at least one oxo group has its carbonyl carbon not in a carboxylic acid.
    additional_oxo_found = False
    for match in ketone_matches:
        # match[1] is the carbonyl carbon in the SMARTS.
        if match[1] not in terminal_acid_indices:
            additional_oxo_found = True
            break
    if not additional_oxo_found:
        for match in aldehyde_matches:
            if match[1] not in terminal_acid_indices:
                additional_oxo_found = True
                break
    if not additional_oxo_found:
        return False, "No additional oxo (aldehyde/ketone) group detected outside the acid"

    # ---------------------------------------------
    # (C) Check that the molecule is “fatty” via its acyclic carbon chain.
    # We require a minimum number of carbon atoms and a long contiguous acyclic carbon chain.
    total_carbons = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:
            total_carbons += 1
    if total_carbons < 12:
        return False, "Too few carbon atoms to be considered a fatty acid"

    # Build a graph of acyclic carbon atoms (carbons not in rings).
    # Only include carbons that are not in any ring.
    acyclic_carbons = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6 and not atom.IsInRing()]
    if not acyclic_carbons:
        return False, "No acyclic carbon chain detected; expected a long aliphatic chain"

    # Build an adjacency list among these acyclic carbon atoms.
    acyclic_adj = {idx: set() for idx in acyclic_carbons}
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtomIdx()
        a2 = bond.GetEndAtomIdx()
        # Only consider bonds between acyclic carbons.
        if a1 in acyclic_adj and a2 in acyclic_adj:
            acyclic_adj[a1].add(a2)
            acyclic_adj[a2].add(a1)

    # Find connected components among acyclic carbons and compute the "chain length" (diameter) for each component.
    visited = set()
    longest_chain = 0
    for node in acyclic_adj:
        if node in visited:
            continue
        # Get the connected component using DFS.
        component = set()
        stack = [node]
        while stack:
            current = stack.pop()
            if current not in component:
                component.add(current)
                for neighbor in acyclic_adj[current]:
                    if neighbor not in component:
                        stack.append(neighbor)
        visited.update(component)
        # For the (acyclic) connected component, compute the diameter (longest shortest-path distance).
        # Use two passes of BFS (first pass: pick any node; second pass: from farthest node, measure distance).
        def bfs(start, nodes):
            dist = {start: 0}
            queue = [start]
            while queue:
                curr = queue.pop(0)
                for nn in acyclic_adj[curr]:
                    if nn in nodes and nn not in dist:
                        dist[nn] = dist[curr] + 1
                        queue.append(nn)
            return dist
        arbitrary = next(iter(component))
        first_dist = bfs(arbitrary, component)
        farthest = max(first_dist, key=first_dist.get)
        second_dist = bfs(farthest, component)
        comp_diameter = max(second_dist.values())
        if comp_diameter > longest_chain:
            longest_chain = comp_diameter

    # For a fatty acid, expect a long chain. Here we require a minimum chain length of 8 bonds (9 carbons).
    if longest_chain < 8:
        return False, "Longest contiguous acyclic carbon chain is too short to be a fatty acid"
    
    # ---------------------------------------------
    # Passed all tests.
    return True, ("Contains a terminal carboxylic acid group, an additional oxo group "
                  "outside that acid, and a sufficiently long acyclic carbon chain typical of fatty acids")

# Example usage:
if __name__ == "__main__":
    # One test example: (6E,8Z)-5-oxooctadecadienoic acid
    test_smiles = "O=C(CCCC(O)=O)/C=C/C=C\\CCCCCCCCC"
    result, reason = is_oxo_fatty_acid(test_smiles)
    print("Test molecule classification:", result)
    print("Reason:", reason)