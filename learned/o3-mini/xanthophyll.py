"""
Classifies: CHEBI:27325 xanthophyll
"""
"""
Classifies: Xanthophyll
Definition: Xanthophylls are oxygenated carotenoids having a long conjugated polyene chain,
an extensive carbon skeleton (typically ≥40 carbons) and at least one oxygen atom.
Many xanthophylls have cyclic endgroups (and then a conjugated chain threshold of ≥4 double bonds),
while acyclic xanthophylls tend to have even longer conjugated chains (≥5 consecutive conjugated double bonds).
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def longest_conjugated_chain(mol):
    """
    Finds the longest chain (i.e. path) of consecutive bonds in the molecule that are
    both double bonds and flagged as conjugated.
    Returns the maximum number of consecutive conjugated double bonds.
    
    We do this by first collecting all bonds that satisfy the condition,
    then building a graph where each bond is a node and two bonds are connected if they share an atom.
    A DFS (depth-first search) is then used to find the longest simple path.
    """
    # collect indices for bonds that are double and conjugated
    conc_dbl_bond_indices = []
    for bond in mol.GetBonds():
        if bond.GetBondTypeAsDouble() == 2.0 and bond.GetIsConjugated():
            conc_dbl_bond_indices.append(bond.GetIdx())
    
    if not conc_dbl_bond_indices:
        return 0

    # Build a mapping from bond index to the two atom indices it connects.
    bond_atoms = {}
    for bond in mol.GetBonds():
        idx = bond.GetIdx()
        if idx in conc_dbl_bond_indices:
            bond_atoms[idx] = (bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
    
    # Build a graph among these bonds: if two bonds share an atom, add an edge (undirected).
    graph = {bond_idx: [] for bond_idx in conc_dbl_bond_indices}
    for i in conc_dbl_bond_indices:
        a1, a2 = bond_atoms[i]
        for j in conc_dbl_bond_indices:
            if i == j:
                continue
            b1, b2 = bond_atoms[j]
            # if any atom is shared then add connection
            if a1 in (b1, b2) or a2 in (b1, b2):
                graph[i].append(j)
    
    # Now perform a DFS for each bond to compute the longest path length in this bond graph.
    visited_global = {}
    def dfs(bond_idx, visited):
        max_length = 1  # count current bond as length 1
        for neighbor in graph[bond_idx]:
            if neighbor not in visited:
                visited.add(neighbor)
                length = 1 + dfs(neighbor, visited)
                if length > max_length:
                    max_length = length
                visited.remove(neighbor)
        return max_length

    max_chain = 0
    for bond_idx in conc_dbl_bond_indices:
        visited = set([bond_idx])
        chain_len = dfs(bond_idx, visited)
        if chain_len > max_chain:
            max_chain = chain_len

    return max_chain

def is_xanthophyll(smiles: str):
    """
    Determines if a molecule is a xanthophyll based on its SMILES string.
    
    A xanthophyll is defined here as an oxygenated carotenoid. The molecule must satisfy:
      - At least one oxygen atom.
      - A carbon scaffold with at least 40 carbon atoms.
      - A minimum molecular weight of 400 Da.
      - A long conjugated polyene chain:
          * If a ring is present, the longest consecutive chain of conjugated double bonds must be ≥4.
          * If acyclic, then it must be ≥5.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is likely a xanthophyll, False otherwise.
        str: A reason for the classification decision.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count oxygen atoms
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 1:
        return False, "No oxygen atoms found; not an oxygenated carotene"

    # Count carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 40:
        return False, f"Too few carbon atoms ({c_count}); expected at least 40 for a carotenoid backbone"

    # Check molecular weight
    mw = rdMolDescriptors.CalcExactMolWt(mol)
    if mw < 400:
        return False, f"Molecular weight too low ({mw:.1f} Da); expected ≥400 Da for a xanthophyll"

    # Check for rings (many xanthophylls have cyclic endgroups, but not all)
    rings = mol.GetRingInfo().AtomRings()
    has_ring = len(rings) > 0

    # Compute the longest consecutive conjugated chain of double bonds
    longest_chain = longest_conjugated_chain(mol)
    # Set the threshold: if ring exists, require at least 4 conjugated double bonds;
    # otherwise for acyclic structures, require at least 5.
    required_chain = 4 if has_ring else 5
    if longest_chain < required_chain:
        chain_req_text = (f"at least {required_chain} consecutive conjugated double bonds "
                          f"(found {longest_chain})")
        return False, f"No long conjugated polyene chain found; expected {chain_req_text}"

    ring_info = "with at least one ring" if has_ring else "without ring structures (acyclic)"
    reason = (f"Contains {ring_info}, a longest conjugated polyene chain of {longest_chain} double bonds, "
              f"sufficient carbon scaffold ({c_count} C atoms) and oxygen functionalities (O count: {o_count}), "
              f"with molecular weight {mw:.1f} Da consistent with xanthophylls")
    return True, reason

# Example usage (for testing); remove or comment out before production deployment:
if __name__ == '__main__':
    test_examples = [
        # Cyclic example: (5R,5'R,6S,8'R)-Luteochrome
        "O1C2(C1(CCCC2(C)C)\\C=C/C(/C)=C\\C=C\\C(\\C)=C\\C=C/C=C(/C=C/C=C(\\C3OC4(C(C(CCC4)(C)C)=C3)C)/C)\\C",
        # Acyclic candidate: spheroidene
        "COC(C)(C)C\\C=C\\C(C)=C\\C=C\\C(C)=C\\C=C\\C(C)=C\\C=C\\C=C(C)\\C=C\\C=C(/C)CC\\C=C(/C)CCC=C(C)C",
    ]
    for sm in test_examples:
        result, reason = is_xanthophyll(sm)
        print("SMILES:", sm)
        print("Result:", result, "|", reason)
        print("-" * 80)