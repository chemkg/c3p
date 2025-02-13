"""
Classifies: CHEBI:10036 wax ester
"""
"""
Classifies: Wax Ester
Definition: A fatty acid ester resulting from the condensation of the carboxy group 
of a fatty acid with the alcoholic hydroxy group of a fatty alcohol.
This module attempts to ensure that (a) exactly one ester group is present, 
(b) the two fragments connected by that group are “simple” (acyclic and solely composed of carbon atoms),
and (c) their sizes (number of carbon atoms) are within a reasonable range.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_wax_ester(smiles: str):
    """
    Determines if a molecule is a simple wax ester based on its SMILES string.
    A wax ester is defined as a fatty acid ester (one ester linkage connecting a fatty acid and a fatty alcohol).
    This function uses an ester SMARTS to identify the single ester bond and then breaks the molecule 
    to extract the two fragments. It then verifies that the fragments contain long, acyclic, predominately carbon chains.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule is classified as a wax ester, False otherwise.
        str: Reason for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Use a SMARTS that matches a simple ester group: C(=O)O (the O has no hydrogen)
    ester_smarts = Chem.MolFromSmarts("[CX3](=O)[OX2H0]")
    ester_matches = mol.GetSubstructMatches(ester_smarts)
    if len(ester_matches) == 0:
        return False, "No ester group found that meets the wax ester pattern"
    if len(ester_matches) > 1:
        return False, "More than one ester group found; not a simple wax ester"
    
    # The ester pattern [CX3](=O)[OX2H0] returns two atoms: index0 = acid carbon (carbonyl) and index1 = ester oxygen.
    acid_c_idx, ester_oxygen_idx = ester_matches[0]
    
    # From the ester oxygen, find the neighboring atom that is not the acid carbon.
    oxy_atom = mol.GetAtomWithIdx(ester_oxygen_idx)
    alcohol_neighbor = None
    for nbr in oxy_atom.GetNeighbors():
        if nbr.GetIdx() != acid_c_idx:
            alcohol_neighbor = nbr
            break
    if alcohol_neighbor is None:
        return False, "Could not determine fatty alcohol side of the ester"
    alcohol_c_idx = alcohol_neighbor.GetIdx()
    
    # Check that the ester group itself is not in a ring
    if oxy_atom.IsInRing() or mol.GetAtomWithIdx(acid_c_idx).IsInRing():
        return False, "Ester group is involved in a ring and thus not a simple wax ester"
    
    # Helper function:
    # DFS over the molecule restricted only to carbon atoms (atomic number 6) and
    # also require that none of the atoms is in a ring.
    def dfs_carbon(start_idx, forbidden):
        visited = set()
        stack = [start_idx]
        while stack:
            curr = stack.pop()
            if curr in visited:
                continue
            # Only traverse if the atom is carbon and not forbidden.
            atom = mol.GetAtomWithIdx(curr)
            if atom.GetAtomicNum() != 6:
                continue
            # Also if any carbon is in a ring, do not include it (expect wax ester chains to be acyclic).
            if atom.IsInRing():
                continue
            visited.add(curr)
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() not in visited and nbr.GetIdx() not in forbidden:
                    stack.append(nbr.GetIdx())
        return visited

    # For the fatty acid portion, start from the acid carbon.
    # We "block" the ester oxygen so that we don't cross back over.
    acid_fragment = dfs_carbon(acid_c_idx, forbidden={ester_oxygen_idx})
    # For the fatty alcohol portion, start from the alcohol carbon.
    alcohol_fragment = dfs_carbon(alcohol_c_idx, forbidden={ester_oxygen_idx})
    
    MIN_CARBONS = 6  # minimum number of carbon atoms required per chain (a low threshold)
    if len(acid_fragment) < MIN_CARBONS:
        return False, f"Fatty acid portion too short; found {len(acid_fragment)} carbon(s)"
    if len(alcohol_fragment) < MIN_CARBONS:
        return False, f"Fatty alcohol portion too short; found {len(alcohol_fragment)} carbon(s)"
    
    # (Optional) To ensure chain linearity, we compute the longest contiguous carbon chain in each fragment.
    # For each fragment, we build a subgraph of just its carbon atom indices and evaluate the maximum distance between any two atoms (chain diameter).
    def longest_linear_chain(carbon_idxs):
        # Build an adjacency dictionary only among the given indices.
        adj = {idx: [] for idx in carbon_idxs}
        for idx in carbon_idxs:
            atom = mol.GetAtomWithIdx(idx)
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in carbon_idxs:
                    adj[idx].append(nbr.GetIdx())
        # For each node, run a DFS to get the maximum distance.
        def dfs_length(start, current, seen):
            max_len = 0
            for neighbor in adj[current]:
                if neighbor not in seen:
                    new_length = 1 + dfs_length(start, neighbor, seen | {neighbor})
                    if new_length > max_len:
                        max_len = new_length
            return max_len
        max_chain = 0
        for idx in carbon_idxs:
            chain_len = dfs_length(idx, idx, {idx})
            if chain_len > max_chain:
                max_chain = chain_len
        # The number of carbons in the longest chain is chain length + 1.
        return max_chain + 1

    acid_chain_length = longest_linear_chain(acid_fragment)
    alcohol_chain_length = longest_linear_chain(alcohol_fragment)
    
    # Check that the chains are long enough typical of fatty chains.
    MIN_LINEAR_CHAIN = 6
    if acid_chain_length < MIN_LINEAR_CHAIN:
        return False, f"Fatty acid chain is not linear enough; longest chain has {acid_chain_length} C(s)"
    if alcohol_chain_length < MIN_LINEAR_CHAIN:
        return False, f"Fatty alcohol chain is not linear enough; longest chain has {alcohol_chain_length} C(s)"
    
    # Return a diagnosis string stating the carbon counts (using the DFS total counts as approximate chain sizes)
    message = ("CORRECT Wax ester detected: fatty acid chain with {} C(s) and fatty alcohol chain with {} C(s) "
               "connected via a single ester linkage.").format(len(acid_fragment), len(alcohol_fragment))
    return True, message

# Example usage:
if __name__ == "__main__":
    # Try one of the valid wax ester examples:
    test_smiles = "O(CCCCCCCC/C=C\\CCCCCC)C(=O)CCCCCCCCCCCCCCCCCCC"  # Palmitoleyl arachidate
    result, message = is_wax_ester(test_smiles)
    print(result, message)