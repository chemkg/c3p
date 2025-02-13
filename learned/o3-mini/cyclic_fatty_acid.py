"""
Classifies: CHEBI:59238 cyclic fatty acid
"""
"""
Classifies: Cyclic Fatty Acid
Definition: Any fatty acid containing anywhere in its structure a ring of atoms.
A fatty acid is here defined as a molecule containing a carboxylic acid group and a largely aliphatic (carbon‐rich) region.
To reduce false positives (e.g. peptides with carboxylic acid groups and rings) we check for a long contiguous aliphatic chain (non‐ring, non‐aromatic carbon atoms).
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def longest_aliphatic_chain(mol):
    """
    Computes the length (number of atoms) of the longest contiguous chain
    composed of aliphatic (sp3) carbon atoms that are not in any ring.
    """
    # First, identify indices of eligible atoms: carbon, non‐aromatic and not in a ring.
    eligible_idxs = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6 and not atom.GetIsAromatic() and not atom.IsInRing():
            # Optionally, we can enforce sp3 hybridization but in many cases the absence of aromaticity is enough.
            if atom.GetHybridization() == Chem.rdchem.HybridizationType.SP3:
                eligible_idxs.append(atom.GetIdx())
    if not eligible_idxs:
        return 0

    # Build an undirected graph: for every eligible atom, list its eligible neighbors.
    graph = {idx: [] for idx in eligible_idxs}
    for idx in eligible_idxs:
        atom = mol.GetAtomWithIdx(idx)
        for neighbor in atom.GetNeighbors():
            n_idx = neighbor.GetIdx()
            if n_idx in graph:  # neighbor is also eligible
                graph[idx].append(n_idx)

    # Given our small graphs, we can use DFS to compute the longest simple path.
    # (Longest path is NP-hard, but our “chain‐graph” here is small.)
    max_chain = 0
    def dfs(current, visited):
        nonlocal max_chain
        visited.add(current)
        current_length = 1  # count the current atom
        branch_max = 0
        for neighbor in graph[current]:
            if neighbor not in visited:
                path_length = dfs(neighbor, visited)
                if path_length > branch_max:
                    branch_max = path_length
        visited.remove(current)
        total = current_length + branch_max
        if total > max_chain:
            max_chain = total
        return total

    # Run DFS starting from every eligible atom.
    for idx in eligible_idxs:
        dfs(idx, set())
    return max_chain

def is_cyclic_fatty_acid(smiles: str):
    """
    Determines if a molecule is a cyclic fatty acid.
    A cyclic fatty acid is defined as a fatty acid (a molecule with a carboxylic acid group and a long aliphatic region)
    that also contains at least one ring anywhere in its structure.
    
    In addition to checking for a carboxylic acid group (SMARTS "C(=O)[O;H]") and a minimum number of carbons,
    we require that there be a contiguous aliphatic chain (non‐aromatic and non‐ring carbon atoms) of a given minimum length.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule is classified as a cyclic fatty acid, False otherwise.
        str: Reason for classification.
    """
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for the presence of a carboxylic acid group.
    acid_pattern = Chem.MolFromSmarts("C(=O)[O;H]")
    if not mol.HasSubstructMatch(acid_pattern):
        return False, "No carboxylic acid group detected; not a fatty acid"
    
    # Count total carbon atoms.
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 7:
        return False, f"Too few carbons (found {c_count}) for a fatty acid"
    
    # Check if the molecule contains at least one ring.
    if mol.GetRingInfo().NumRings() < 1:
        return False, "No ring detected in the structure; not a cyclic fatty acid"
    
    # Check for a long contiguous aliphatic chain.
    chain_length = longest_aliphatic_chain(mol)
    # Use a threshold (arbitrarily 6 carbon atoms) to ensure a significant aliphatic region.
    if chain_length < 6:
        return False, f"Longest contiguous aliphatic chain is too short (length {chain_length}); not a fatty acid"
    
    return True, "Contains carboxylic acid group, a ring, and a long aliphatic chain; qualifies as a cyclic fatty acid"

# Example testing (can be removed or commented out in production)
if __name__ == "__main__":
    test_smiles = [
        "O=C(CCC/C=C\\C/C=C\\C/C=C\\[C@@H]([C@H]1[C@H](CCCCC)O1)O)O",  # cyclic fatty acid true positive
        "O=C(N1[C@H](CCC1)C(O)=O)[C@@H](NC(=O)[C@@H](N)CC(O)=O)CCC(=O)N"   # peptide false positive
    ]
    for s in test_smiles:
        decision, reason = is_cyclic_fatty_acid(s)
        print(f"SMILES: {s}\nResult: {decision}\nReason: {reason}\n")