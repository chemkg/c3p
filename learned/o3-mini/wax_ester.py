"""
Classifies: CHEBI:10036 wax ester
"""
"""
Classifies: Wax Ester
Definition: A fatty acid ester resulting from the condensation of the carboxy group 
of a fatty acid with the alcoholic hydroxy group of a fatty alcohol.
This implementation ensures that exactly one ester group is present, 
extracts the two fragments attached to that ester, and checks whether they are long, 
acyclic and primarily composed of carbon atoms.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_wax_ester(smiles: str):
    """
    Determines if a molecule is a simple wax ester based on its SMILES string.
    A wax ester is defined as a fatty acid ester (one ester linkage connecting a fatty acid and a fatty alcohol).
    This function uses an ester SMARTS to identify the single ester bond and then partitions the molecule 
    into two fragments. It then verifies that the fragments are long carbon chains (acyclic and predominately carbon).

    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule is classified as a wax ester, False otherwise.
        str: Reason for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Use a SMARTS that matches a simple ester group: C(=O)O 
    # The pattern "[CX3](=O)[OX2H0]" will match three atoms: 
    # index 0: the acid carbon, index 1: the carbonyl oxygen, index 2: the ester oxygen.
    ester_smarts = Chem.MolFromSmarts("[CX3](=O)[OX2H0]")
    ester_matches = mol.GetSubstructMatches(ester_smarts)
    if len(ester_matches) == 0:
        return False, "No ester group found that meets the wax ester pattern"
    if len(ester_matches) > 1:
        return False, "More than one ester group found; not a simple wax ester"
    
    # Unpack the matched atoms: acid carbon is at index 0 and ester oxygen at index 2.
    match = ester_matches[0]
    acid_c_idx = match[0]
    ester_oxygen_idx = match[2]
    
    # From the ester oxygen, determine the fatty alcohol portion.
    oxy_atom = mol.GetAtomWithIdx(ester_oxygen_idx)
    alcohol_neighbor = None
    for nbr in oxy_atom.GetNeighbors():
        if nbr.GetIdx() != acid_c_idx:
            alcohol_neighbor = nbr
            break
    if alcohol_neighbor is None:
        return False, "Could not determine fatty alcohol side of the ester"
    alcohol_c_idx = alcohol_neighbor.GetIdx()
    
    # Check that the ester group is not in a ring (we expect a simple linear ester linkage).
    if oxy_atom.IsInRing() or mol.GetAtomWithIdx(acid_c_idx).IsInRing():
        return False, "Ester group is involved in a ring and thus not a simple wax ester"
    
    # Helper function: perform a depth-first search (DFS) from a starting atom and collect connected carbons.
    # We only traverse through carbon atoms (atomic number 6) that are not in rings.
    def dfs_carbon(start_idx, forbidden):
        visited = set()
        stack = [start_idx]
        while stack:
            curr = stack.pop()
            if curr in visited:
                continue
            atom = mol.GetAtomWithIdx(curr)
            if atom.GetAtomicNum() != 6:
                continue
            if atom.IsInRing():
                continue
            visited.add(curr)
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() not in visited and nbr.GetIdx() not in forbidden:
                    stack.append(nbr.GetIdx())
        return visited

    # For the fatty acid portion, start from the acid carbon.
    # Block the traversal from crossing back over the ester oxygen.
    acid_fragment = dfs_carbon(acid_c_idx, forbidden={ester_oxygen_idx})
    # For the fatty alcohol portion, start from the alcohol carbon.
    alcohol_fragment = dfs_carbon(alcohol_c_idx, forbidden={ester_oxygen_idx})
    
    MIN_CARBONS = 6  # minimum number of carbon atoms required per chain
    if len(acid_fragment) < MIN_CARBONS:
        return False, f"Fatty acid portion too short; found {len(acid_fragment)} carbon(s)"
    if len(alcohol_fragment) < MIN_CARBONS:
        return False, f"Fatty alcohol portion too short; found {len(alcohol_fragment)} carbon(s)"
    
    # (Optional) To ensure chain linearity, compute the longest contiguous chain of carbons.
    def longest_linear_chain(carbon_idxs):
        # Build an adjacency list among the provided indices.
        adj = {idx: [] for idx in carbon_idxs}
        for idx in carbon_idxs:
            atom = mol.GetAtomWithIdx(idx)
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in carbon_idxs:
                    adj[idx].append(nbr.GetIdx())
        # DFS to find maximum chain length.
        def dfs_length(current, seen):
            max_len = 0
            for neighbor in adj[current]:
                if neighbor not in seen:
                    new_length = 1 + dfs_length(neighbor, seen | {neighbor})
                    if new_length > max_len:
                        max_len = new_length
            return max_len
        max_chain = 0
        for idx in carbon_idxs:
            chain_len = dfs_length(idx, {idx})
            if chain_len > max_chain:
                max_chain = chain_len
        # Longest chain length in terms of number of atoms.
        return max_chain + 1

    acid_chain_length = longest_linear_chain(acid_fragment)
    alcohol_chain_length = longest_linear_chain(alcohol_fragment)
    
    MIN_LINEAR_CHAIN = 6
    if acid_chain_length < MIN_LINEAR_CHAIN:
        return False, f"Fatty acid chain is not linear enough; longest chain has {acid_chain_length} C(s)"
    if alcohol_chain_length < MIN_LINEAR_CHAIN:
        return False, f"Fatty alcohol chain is not linear enough; longest chain has {alcohol_chain_length} C(s)"
    
    message = ("CORRECT Wax ester detected: fatty acid chain with {} C(s) and fatty alcohol chain with {} C(s) "
               "connected via a single ester linkage.").format(len(acid_fragment), len(alcohol_fragment))
    return True, message

# Example usage:
if __name__ == "__main__":
    # Test one of the provided wax ester examples:
    test_smiles = "O(CCCCCCCCCCCCCCCCCCCC(C)C)C(CCCCCCC/C=C\\CCCCCCCC)=O"  # 1-O-20-methylhenicosyl oleate
    result, message = is_wax_ester(test_smiles)
    print(result, message)