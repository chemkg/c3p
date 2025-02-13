"""
Classifies: CHEBI:35366 fatty acid
"""
"""
Classifies: Fatty Acids 
Definition: “Fatty acids” here are defined as acyclic aliphatic monocarboxylic acids.
They must be a valid molecule with no rings, contain exactly one terminal carboxyl group 
(a carbonyl carbon with one –OH or [O-] attached that is attached to only one carbon),
have no amide bond(s) (to avoid peptides), and show a predominantly linear hydrocarbon chain.
In addition we require that the longest continuous carbon chain (backbone) constitutes a large 
fraction of the total carbons and that the ratio of extra heteroatoms (aside from the two in the acid)
relative to carbon count is not excessive.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_fatty_acid(smiles: str):
    """
    Determines if a molecule is a fatty acid based on its SMILES string.
    
    The molecule must:
      1. Be a valid molecule.
      2. Be acyclic (no rings).
      3. Contain exactly one terminal carboxylic acid group – a carbon atom
         with a double-bonded oxygen and a singly-bonded hydroxyl (or O-),
         where that acid carbon is attached to exactly one other carbon.
      4. Not contain amide bonds (e.g. C(=O)N).
      5. Have a sufficiently long main carbon chain that accounts for most of the molecule’s carbons.
      6. Not contain an excess of heteroatoms (or disallowed atoms such as phosphorus) relative to a typical fatty acid.
      
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a fatty acid, False otherwise.
        str: A textual explanation for the decision.
    """
    # Parse the molecule from SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Reject if the molecule contains any rings
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule contains ring(s), expected an acyclic fatty acid"
    
    # Define SMARTS for a carboxylic acid group.
    # The pattern matches a carbonyl carbon (C(=O)) with a singly-bound oxygen that is either protonated (–OH) or deprotonated ([O-]).
    carboxyl_smarts = "C(=O)[O;H1,O-]"
    carboxyl_query = Chem.MolFromSmarts(carboxyl_smarts)
    matches = mol.GetSubstructMatches(carboxyl_query)
    
    # Ensure that the carboxyl group is terminal.
    # For each match, check that the carboxyl carbon is attached to exactly one carbon.
    terminal_matches = []
    for match in matches:
        c_idx = match[0]  # the carbonyl carbon from the carboxyl group
        c_atom = mol.GetAtomWithIdx(c_idx)
        # Count how many of its neighbors are carbon atoms.
        carbon_neighbors = [nbr for nbr in c_atom.GetNeighbors() if nbr.GetAtomicNum() == 6]
        if len(carbon_neighbors) == 1:
            terminal_matches.append(match)
    
    if len(terminal_matches) != 1:
        return False, f"Found {len(terminal_matches)} terminal carboxyl group(s), expected exactly 1"
    
    # Reject molecules with amide bonds which are common in peptides
    amide_query = Chem.MolFromSmarts("C(=O)N")
    if mol.HasSubstructMatch(amide_query):
        return False, "Molecule contains amide bond(s), likely not a free fatty acid"
    
    # Reject if disallowed atoms (like phosphorus) are present.
    p_query = Chem.MolFromSmarts("[#15]")
    if mol.HasSubstructMatch(p_query):
        return False, "Molecule contains phosphorus, likely not a free fatty acid"
    
    # Count total number of carbon atoms.
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 4:
        return False, f"Too few carbon atoms ({carbon_count}) to be a fatty acid"
    
    # Build a graph over only the carbon atoms.
    carbon_indices = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    graph = {idx: [] for idx in carbon_indices}
    for bond in mol.GetBonds():
        a = bond.GetBeginAtomIdx()
        b = bond.GetEndAtomIdx()
        if a in graph and b in graph:
            graph[a].append(b)
            graph[b].append(a)
    
    # Identify the carbon neighbor of the terminal carboxyl group.
    # (In our terminal match, match[0] is the carboxyl carbon.)
    carboxyl_carbon_idx = terminal_matches[0][0]
    carboxyl_atom = mol.GetAtomWithIdx(carboxyl_carbon_idx)
    neighbors = [nbr.GetIdx() for nbr in carboxyl_atom.GetNeighbors() if nbr.GetAtomicNum() == 6]
    if not neighbors:
        return False, "Terminal carboxyl group has no attached carbon"
    chain_start = neighbors[0]
    
    # Compute the longest chain (i.e. the diameter) of the carbon subgraph in the connected component that contains chain_start.
    # Because the molecule is acyclic, the carbon–only subgraph is a tree.
    from collections import deque
    def bfs(start):
        visited = {start}
        queue = deque([(start, 0)])
        farthest = (start, 0)
        while queue:
            node, dist = queue.popleft()
            if dist > farthest[1]:
                farthest = (node, dist)
            for nb in graph[node]:
                if nb not in visited:
                    visited.add(nb)
                    queue.append((nb, dist + 1))
        return farthest, visited

    # Get the connected component (set of carbon indices) that includes chain_start.
    _, comp = bfs(chain_start)
    any_node = next(iter(comp))
    (node_A, _), _ = bfs(any_node)
    (node_B, diameter), _ = bfs(node_A)
    longest_chain_length = diameter + 1  # number of carbons in the longest path
    
    # Require that the longest carbon chain (backbone) constitutes most of the molecule.
    ratio = longest_chain_length / carbon_count
    # We require at least 70% of all carbons to be in one contiguous chain.
    if ratio < 0.7:
        return False, f"Longest carbon chain length ({longest_chain_length}) is too short relative to total carbons ({carbon_count})"
    
    # Heuristic on heteroatoms:
    # Count heteroatoms (all atoms except C and H)
    hetero_total = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() not in (1, 6))
    # The carboxyl group brings 2 oxygens that are expected.
    extra_hetero = hetero_total - 2  
    # Allow roughly one extra heteroatom per 5 carbons; adjust threshold as needed.
    if extra_hetero > (carbon_count // 5):
        return False, f"High heteroatom content (extra={extra_hetero} for {carbon_count} carbons), not typical for a fatty acid"
    
    return True, f"Molecule is an acyclic aliphatic carboxylic acid with {carbon_count} carbon(s) and a main chain of {longest_chain_length} carbon(s)"

# Example usage:
if __name__ == "__main__":
    # Test with butyric acid SMILES (should be True)
    test_smiles = "CCCC(O)=O"
    result, reason = is_fatty_acid(test_smiles)
    print(result, reason)