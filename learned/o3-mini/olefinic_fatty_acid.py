"""
Classifies: CHEBI:53339 olefinic fatty acid
"""
"""
Classifies: Olefinic Fatty Acid
Definition: Any fatty acid containing at least one C=C double bond.
A valid olefinic fatty acid (for this classifier) is defined here as a long‐chain (>=12 carbon atoms),
acyclic molecule that:
  - Contains exactly one terminal carboxylic acid group (the acid C is connected to only one carbon).
  - Contains at least one C=C double bond.
  - Has a high carbon fraction among non‐hydrogen atoms.
  - Has a nearly linear (unbranched) alkyl chain: the longest continuous carbon chain accounts for
    at least 90% of all carbon atoms.
Additionally, molecules containing phosphorus are rejected (e.g. as in complex lipids).
Note: This heuristic may fail for edge cases.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from collections import deque

def is_olefinic_fatty_acid(smiles: str):
    """
    Determines if a molecule is an olefinic fatty acid based on its SMILES string.
    Olefinic fatty acids are defined as acyclic long‐chain carboxylic acids
    with a terminal carboxylic acid group (the acid carbon being bonded to only one carbon),
    containing at least one C=C double bond, and made mostly of carbon atoms. In addition,
    the continuous carbon chain (calculated as the longest path in the carbon–carbon graph)
    should account for nearly all of the carbon atoms (>=90%), suggesting a simple unadorned chain.

    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule meets the olefinic fatty acid criteria, False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Reject if the molecule contains any rings.
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule contains rings; not a simple acyclic fatty acid"
    
    # Reject if the molecule contains phosphorus (typically in phospholipids).
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 15:
            return False, "Molecule contains phosphorus; likely not a simple fatty acid"
    
    # Look for a carboxylic acid group.
    # SMARTS: C(=O)[O;H1,-] matches a carbonyl C bonded to an -OH (or negative O).
    acid_smarts = "C(=O)[O;H1,-]"
    acid_pattern = Chem.MolFromSmarts(acid_smarts)
    acid_matches = mol.GetSubstructMatches(acid_pattern)
    if not acid_matches:
        return False, "No carboxylic acid group found; not a fatty acid"
    
    # Require exactly one carboxylic acid group.
    if len(acid_matches) != 1:
        return False, f"Found {len(acid_matches)} carboxylic acid groups; expected exactly one"
        
    # Check that the acid is terminal.
    # The acid carbon should have exactly one carbon neighbor.
    terminal_acid_found = False
    acid_match = acid_matches[0]
    acid_c_idx = acid_match[0]  # first atom in the match is the acid carbon.
    acid_atom = mol.GetAtomWithIdx(acid_c_idx)
    # Count carbon (atomic number 6) neighbors.
    carbon_neighbors = [nbr for nbr in acid_atom.GetNeighbors() if nbr.GetAtomicNum() == 6]
    if len(carbon_neighbors) == 1:
        terminal_acid_found = True
    if not terminal_acid_found:
        return False, "Carboxylic acid group is not terminal; not a typical fatty acid"
    
    # Count total number of carbon atoms.
    carbon_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    n_carbons = len(carbon_atoms)
    if n_carbons < 12:
        return False, f"Too few carbon atoms ({n_carbons} found, need at least 12) to be a long-chain fatty acid"
    
    # Check for the presence of a carbon–carbon double bond (C=C).
    olefin_pattern = Chem.MolFromSmarts("C=C")
    if not mol.HasSubstructMatch(olefin_pattern):
        return False, "No C=C double bonds found; not an olefinic fatty acid"
    
    # Heuristic: the majority of heavy (non-hydrogen) atoms should be carbon.
    non_H_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() != 1]
    n_nonH = len(non_H_atoms)
    if n_nonH == 0 or (n_carbons / n_nonH) < 0.7:
        return False, "The molecule does not appear to be a simple fatty acyl chain (low fraction of carbon atoms)"
    
    # Now check the linearity of the molecule by examining the longest continuous path among carbons.
    # Build a graph (as a dictionary) from the indices of carbon atoms.
    carbon_indices = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    carbon_graph = {idx: [] for idx in carbon_indices}
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        if a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 6:
            idx1 = a1.GetIdx()
            idx2 = a2.GetIdx()
            if idx1 in carbon_graph and idx2 in carbon_graph:
                carbon_graph[idx1].append(idx2)
                carbon_graph[idx2].append(idx1)
                
    # If for some reason we have no carbon–carbon bonds, reject.
    if not carbon_graph:
        return False, "No carbon–carbon bonds found; not a fatty acid"
    
    # Helper: perform BFS from a start node on the carbon graph to get distances.
    def bfs_longest(start, graph):
        visited = {start}
        queue = deque([(start, 0)])
        farthest_node = start
        max_dist = 0
        while queue:
            node, dist = queue.popleft()
            if dist > max_dist:
                max_dist = dist
                farthest_node = node
            for nbr in graph[node]:
                if nbr not in visited:
                    visited.add(nbr)
                    queue.append((nbr, dist+1))
        return farthest_node, max_dist

    # Pick an arbitrary carbon as starting point.
    start_node = next(iter(carbon_graph))
    far_node, _ = bfs_longest(start_node, carbon_graph)
    # Now from far_node, get the maximum distance.
    _, longest_path_length = bfs_longest(far_node, carbon_graph)
    # Since distance counts edges, the number of atoms in the longest chain is (longest_path_length + 1)
    longest_chain_atoms = longest_path_length + 1
    
    # For a simple fatty acid, nearly all carbon atoms should be in the main chain.
    # If the ratio of (longest chain carbons / total carbons) is below 0.9, the chain is too branched.
    if (longest_chain_atoms / n_carbons) < 0.9:
        return False, f"Longest continuous carbon chain accounts for only {longest_chain_atoms} of {n_carbons} carbons; likely branched or decorated"
    
    return True, ("Contains a terminal carboxylic acid group, is acyclic and long‐chained (>=12 C atoms) with a nearly linear "
                  f"alkyl chain ({longest_chain_atoms}/{n_carbons} carbons in the main chain) and has at least one C=C double bond, "
                  "classifying it as an olefinic fatty acid.")

# Example usage:
if __name__ == '__main__':
    # An example known olefinic fatty acid
    test_smiles = "CCCCCC\\C=C/CC\\C=C/CCCC(O)=O"  # a simple acyclic, unsaturated fatty acid
    result, reason = is_olefinic_fatty_acid(test_smiles)
    print(result, reason)