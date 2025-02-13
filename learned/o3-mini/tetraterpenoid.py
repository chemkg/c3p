"""
Classifies: CHEBI:26935 tetraterpenoid
"""
"""
Classifies: Tetraterpenoid
Defined as any terpenoid derived from a tetraterpene (typically with a C40 core,
possibly rearranged or slightly modified). They are generally characterized by an extended
conjugated polyene chain.
This implementation improves the previous version by memoizing the DFS search for the longest
conjugated carbon chain, thereby avoiding a timeout on complex molecules.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from functools import lru_cache

def is_tetraterpenoid(smiles: str):
    """
    Determines if a molecule is a tetraterpenoid based on its SMILES string.
    Improved classification consists of:
      1. Check that the molecular weight is within an expected range (300-750 Da).
      2. Build a carbon-backbone graph including only carbon atoms connected via conjugated bonds.
      3. Find the longest simple path in that graph using a memoized DFS with a bitmask for the visited set.
         The path is taken as a proxy for the tetraterpene-derived conjugated backbone.
      4. Verify that the backbone has a length between 30 and 50 atoms.
      5. Verify that along the backbone there are at least 7 non‐aromatic C=C bonds.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a tetraterpenoid, False otherwise.
        str: Explanation for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # 1. Check overall molecular weight.
    mw = rdMolDescriptors.CalcExactMolWt(mol)
    if mw < 300 or mw > 750:
        return False, f"Molecular weight ({mw:.1f}) is not in the expected range (300-750 Da) for tetraterpenoids."
    
    # 2. Build a conjugated carbon backbone graph.
    # Include only carbon atoms (atomic number 6) that are connected by conjugated bonds.
    backbone_graph = {}
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:
            idx = atom.GetIdx()
            backbone_graph[idx] = []  # initialize neighbor list

    for bond in mol.GetBonds():
        # Only use bond if it is conjugated.
        if bond.GetIsConjugated():
            a1 = bond.GetBeginAtom()
            a2 = bond.GetEndAtom()
            if a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 6:
                i1, i2 = a1.GetIdx(), a2.GetIdx()
                if i1 in backbone_graph and i2 in backbone_graph:
                    backbone_graph[i1].append(i2)
                    backbone_graph[i2].append(i1)
    
    if not backbone_graph:
        return False, "No suitable conjugated carbon substructure found."
    
    # Remap node indices to bit positions for memoization.
    # Create list of nodes in the graph.
    nodes = list(backbone_graph.keys())
    node_to_bit = {node: i for i, node in enumerate(nodes)}
    bit_to_node = {i: node for node, i in node_to_bit.items()}
    n_nodes = len(nodes)
    
    # 3. Use DFS with memoization over bitmask-encoded visited set.
    # Returns a tuple (max_length, best_path) where best_path is a tuple
    # of bit positions representing the longest chain starting from the current node.
    @lru_cache(maxsize=None)
    def dfs(current_bit, visited_mask):
        best_length = 1
        best_path = (current_bit,)
        current_node = bit_to_node[current_bit]
        for neighbor in backbone_graph[current_node]:
            neighbor_bit = node_to_bit[neighbor]
            if not (visited_mask & (1 << neighbor_bit)):
                new_mask = visited_mask | (1 << neighbor_bit)
                candidate_length, candidate_path = dfs(neighbor_bit, new_mask)
                if 1 + candidate_length > best_length:
                    best_length = 1 + candidate_length
                    best_path = (current_bit,) + candidate_path
        return best_length, best_path
    
    # 4. Find the overall longest path in the graph.
    overall_best_length = 0
    overall_best_path = ()
    for bit in range(n_nodes):
        length, path = dfs(bit, 1 << bit)
        if length > overall_best_length:
            overall_best_length = length
            overall_best_path = path
    
    # Check if the longest conjugated carbon chain is in the expected range.
    if overall_best_length < 30 or overall_best_length > 50:
        return (False, 
                f"Longest conjugated carbon chain length ({overall_best_length}) is not in the expected range (30-50) for tetraterpenoids.")
    
    # 5. Count the number of non‐aromatic double bonds along the recovered backbone.
    # Convert best_path (bit indices) back to original atom indices.
    chain_atom_indices = [bit_to_node[b] for b in overall_best_path]
    double_bond_count = 0
    # Walk along consecutive atoms in the chain.
    for i in range(len(chain_atom_indices) - 1):
        bond = mol.GetBondBetweenAtoms(chain_atom_indices[i], chain_atom_indices[i+1])
        if bond is not None and bond.GetBondType() == Chem.BondType.DOUBLE:
            # Exclude aromatic bonds.
            if not (bond.GetBeginAtom().GetIsAromatic() or bond.GetEndAtom().GetIsAromatic()):
                double_bond_count += 1

    if double_bond_count < 7:
        return False, (f"Found only {double_bond_count} non-aromatic C=C bonds along the backbone; "
                       "expected at least 7 for an extended conjugated system.")
    
    return True, (f"Longest conjugated chain contains {overall_best_length} carbons with "
                  f"{double_bond_count} non-aromatic C=C bonds, and molecular weight ({mw:.1f}) "
                  "is in the expected range for a tetraterpenoid.")

# Example usage:
# Uncomment below lines to test one of the provided SMILES examples.
# test_smiles = "COC(C)(C)C\\C=C\\C(C)=C\\C=C\\C(C)=C\\C=C\\C(C)=C\\C=C\\C=C(C)\\C=C\\C=C(/C)CC\\C=C(/C)CCC=C(C)C"
# result, reason = is_tetraterpenoid(test_smiles)
# print(result, reason)