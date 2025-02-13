"""
Classifies: CHEBI:59238 cyclic fatty acid
"""
"""
Classifies: Cyclic Fatty Acid
Definition: Any fatty acid containing anywhere in its structure a ring of atoms.
A fatty acid is defined here as a molecule containing a carboxylic acid group and a largely aliphatic (carbon‐rich) region.
To further reduce false positives such as peptides with carboxylic acids and complex ring systems, we:
  (a) require a long contiguous chain of aliphatic (sp3 non‐aromatic, non‐ring) carbons (minimum length of 5 atoms), and 
  (b) require that this chain account for at least 30% of all carbon atoms in the molecule.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def longest_aliphatic_chain(mol):
    """
    Computes the length (number of atoms) of the longest contiguous chain
    composed only of aliphatic (sp3) carbon atoms that are not aromatic and not in any ring.
    """
    # Identify indices of eligible atoms: carbon atoms that are non‐aromatic,
    # in sp3 hybridization, and not in any ring.
    eligible_idxs = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6 and not atom.GetIsAromatic() and not atom.IsInRing():
            if atom.GetHybridization() == Chem.rdchem.HybridizationType.SP3:
                eligible_idxs.append(atom.GetIdx())
    if not eligible_idxs:
        return 0

    # Build a graph: for every eligible atom, list its eligible neighbors.
    graph = {idx: [] for idx in eligible_idxs}
    for idx in eligible_idxs:
        atom = mol.GetAtomWithIdx(idx)
        for neighbor in atom.GetNeighbors():
            n_idx = neighbor.GetIdx()
            if n_idx in graph:
                graph[idx].append(n_idx)

    max_chain = 0
    # Use DFS (depth first search) to compute the longest simple path in the graph.
    def dfs(current, visited):
        nonlocal max_chain
        visited.add(current)
        current_length = 1  # count current atom
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

    # Run DFS from every eligible atom.
    for idx in eligible_idxs:
        dfs(idx, set())
    return max_chain

def is_cyclic_fatty_acid(smiles: str):
    """
    Determines if a molecule is a cyclic fatty acid.
    A cyclic fatty acid is defined as a fatty acid (a molecule with a carboxylic acid group and a long aliphatic region)
    that also contains at least one ring anywhere in its structure.

    Beyond detecting a carboxylic acid group and a ring, we further require:
      - that there is a contiguous chain of eligible aliphatic carbons (sp3, non‐aromatic, non‐ring) of at least 5 atoms,
      - and that this chain accounts for at least 30% of the total carbon atoms – to avoid false positives from peptide-like molecules.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is classified as a cyclic fatty acid, False otherwise.
        str: Reason for the classification.
    """
    # Convert SMILES into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for the presence of a carboxylic acid group.
    acid_pattern = Chem.MolFromSmarts("C(=O)[O;H]")
    if not mol.HasSubstructMatch(acid_pattern):
        return False, "No carboxylic acid group detected; not a fatty acid"
    
    # Count total carbon atoms.
    c_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    c_count = len(c_atoms)
    if c_count < 7:
        return False, f"Too few carbons (found {c_count}) for a fatty acid"
    
    # Check that at least one ring is present.
    if mol.GetRingInfo().NumRings() < 1:
        return False, "No ring detected in the structure; not a cyclic fatty acid"
    
    # Compute the longest contiguous aliphatic chain (non‐aromatic, sp3, non‐ring carbons).
    chain_length = longest_aliphatic_chain(mol)
    # Require a minimum chain length of 5 atoms.
    if chain_length < 5:
        return False, f"Longest contiguous aliphatic chain is too short (length {chain_length}); not a fatty acid"
    
    # Also require that this chain represents at least 30% of the total carbon atoms.
    aliphatic_fraction = chain_length / c_count
    if aliphatic_fraction < 0.30:
        return False, (f"Longest aliphatic chain accounts for only {aliphatic_fraction:.2f} "
                       f"of total carbons ({chain_length} of {c_count}); not enough aliphatic region")
                       
    return True, ("Contains carboxylic acid group, a ring, and a sufficiently long and dominant aliphatic chain; "
                  "qualifies as a cyclic fatty acid")

# Example testing below (this section can be removed or commented out in production)
if __name__ == "__main__":
    # Test examples (the provided structures from the outcomes)
    test_examples = [
        # True positives – expected to be classified as cyclic fatty acids.
        "O1C(CCCCCCC(O)=O)=C(C(=C1CCC)C)C",  # 3,4-Dimethyl-5-propyl-2-furanheptanoic acid
        "[C@H]1(CCCCCCCC(=O)O)[C@@H](CCCCCCCC)O1",  # (9S,10R)-epoxyoctadecanoic acid
        "O1C(CCCCCCCC(O)=O)=CC=C1CCCCC",  # 5-Pentyl-2-furanoctanoic acid
        "OC(=O)CCCCCCCCCC[C@@H]1CCC=C1",  # (R)-hydnocarpic acid
        "CCCCCCCCC1=C(CCCCCCC(O)=O)C1",  # malvalic acid
        "OC(=O)CCCCCCCCCCCC[C@H]1CCC=C1",  # (S)-chaulmoogric acid
        # False negatives – cyclic fatty acids that may have been previously missed
        "O=C(CCC/C=C\\C/C=C\\C/C=C\\[C@@H]([C@H]1[C@H](CCCCC)O1)O)O",  # (13S)-hydroxy-(14S,15S)-epoxy-(5Z,8Z,11Z)-icosatrienoic acid
        "C(CCC/C=C\\C[C@@H]1/C(/O1)=C/C=C\\C/C=C\\CCCCC)(=O)O",  # (5Z,8R,9Z,11Z,14Z)-8,9-epoxyicosatetraenoic acid
    ]
    for s in test_examples:
        decision, reason = is_cyclic_fatty_acid(s)
        print(f"SMILES: {s}\nResult: {decision}\nReason: {reason}\n")