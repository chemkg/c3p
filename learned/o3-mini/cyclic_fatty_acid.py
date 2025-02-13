"""
Classifies: CHEBI:59238 cyclic fatty acid
"""
"""
Classifies: Cyclic Fatty Acid
Definition: Any fatty acid containing anywhere in its structure a ring of atoms.
A fatty acid (here) is a molecule that has:
  - A carboxylic acid group ("C(=O)[O;H]")
  - A long, largely aliphatic (carbon‐rich) region. To ensure this we require:
      • the molecule has at least 12 carbon atoms overall,
      • a contiguous chain (of sp3, non‐aromatic, non‐ring) of at least 6 carbon atoms,
      • and that this chain represents at least 25% of the total carbon atoms.
  - At least one ring (anywhere in the structure).
  - Also, to avoid misclassifying peptides or other non–fatty acid molecules, we require that the overall fraction
    of carbons (versus all heavy atoms) is at least 50%.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def longest_aliphatic_chain(mol):
    """
    Computes the length (number of atoms) of the longest contiguous chain composed only of aliphatic (sp3) carbon atoms
    that are non‐aromatic and not in any ring.
    """
    # Get indices of eligible atoms: carbon atoms that are sp3, non‐aromatic and not in any ring.
    eligible_idxs = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6 and not atom.GetIsAromatic() and not atom.IsInRing():
            if atom.GetHybridization() == Chem.rdchem.HybridizationType.SP3:
                eligible_idxs.append(atom.GetIdx())
    if not eligible_idxs:
        return 0

    # Build a graph: each eligible atom and its eligible neighbors.
    graph = {idx: [] for idx in eligible_idxs}
    for idx in eligible_idxs:
        atom = mol.GetAtomWithIdx(idx)
        for neighbor in atom.GetNeighbors():
            n_idx = neighbor.GetIdx()
            if n_idx in graph:
                graph[idx].append(n_idx)

    max_chain = 0
    # Recursive DFS to find the longest simple path.
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

    # Perform DFS from all eligible atoms.
    for idx in eligible_idxs:
        dfs(idx, set())
    return max_chain

def is_cyclic_fatty_acid(smiles: str):
    """
    Determines if a molecule (given its SMILES string) is a cyclic fatty acid.
    A cyclic fatty acid is defined here as a molecule that has:
      • at least one carboxylic acid group (pattern: C(=O)[O;H]),
      • at least one ring anywhere in its structure,
      • an overall long aliphatic region: the molecule must contain at least 12 carbon atoms,
          and a contiguous chain (of sp3, non‐aromatic, non‐ring carbons) of at least 6 atoms that accounts
          for at least 25% of all carbon atoms in the molecule;
      • and the molecule as a whole must be carbon‐rich (at least 50% of its heavy atoms are carbons).
    
    Args:
        smiles (str): SMILES string representing the molecule.
    
    Returns:
        bool: True if the molecule qualifies as a cyclic fatty acid, False otherwise.
        str: Explanation for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group ("C(=O)[O;H]")
    acid_pattern = Chem.MolFromSmarts("C(=O)[O;H]")
    if not mol.HasSubstructMatch(acid_pattern):
        return False, "No carboxylic acid group detected; not a fatty acid"

    # Count total carbon atoms.
    c_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    c_count = len(c_atoms)
    if c_count < 12:
        return False, f"Too few carbons (found {c_count}); fatty acids typically have at least 12 carbons"

    # Check that at least one ring is present in the structure.
    if mol.GetRingInfo().NumRings() < 1:
        return False, "No ring detected; not a cyclic fatty acid"

    # Compute the longest contiguous aliphatic chain.
    chain_length = longest_aliphatic_chain(mol)
    if chain_length < 6:
        return False, f"Longest contiguous aliphatic chain is too short (length {chain_length}); not a fatty acid"
    
    # Require that this chain accounts for at least 25% of total carbons.
    aliphatic_fraction = chain_length / c_count
    if aliphatic_fraction < 0.25:
        return False, (f"Longest aliphatic chain accounts for only {aliphatic_fraction:.2f} of total carbons "
                       f"({chain_length} of {c_count}); not enough aliphatic region")

    # Check overall carbon richness: at least 50% of heavy atoms (atoms with atomic number >1) should be carbon.
    heavy_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1]
    if heavy_atoms:
        carbon_fraction = c_count / len(heavy_atoms)
        if carbon_fraction < 0.5:
            return False, (f"Only {carbon_fraction:.2f} of heavy atoms are carbons; molecule is not sufficiently carbon-rich")
    else:
        return False, "No heavy atoms found; invalid molecule"
    
    return True, ("Contains a carboxylic acid group, a ring, and a sufficiently long aliphatic region "
                  "with a dominant carbon content; qualifies as a cyclic fatty acid")

# The following test code is for demonstration (can be removed or commented out in production).
if __name__ == "__main__":
    test_smiles_list = [
        # True positives (expected to be cyclic fatty acids)
        "O1C(CCCCCCC(O)=O)=C(C(=C1CCC)C)C",  # 3,4-Dimethyl-5-propyl-2-furanheptanoic acid
        "[C@H]1(CCCCCCCC(=O)O)[C@@H](CCCCCCCC)O1",  # (9S,10R)-epoxyoctadecanoic acid
        "O1C(CCCCCCCC(O)=O)=CC=C1CCCCC",  # 5-Pentyl-2-furanoctanoic acid
        "OC(=O)CCCCCCCCCC[C@@H]1CCC=C1",  # (R)-hydnocarpic acid
        "CCCCCCCCC1=C(CCCCCCC(O)=O)C1",  # malvalic acid
        "OC(=O)CCCCCCCCCCCC[C@H]1CCC=C1",  # (S)-chaulmoogric acid
        # False negatives (should be cyclic fatty acids but previously were missed)
        "O=C(CCC/C=C\\C/C=C\\C/C=C\\[C@@H]([C@H]1[C@H](CCCCC)O1)O)O",  # (13S)-hydroxy-(14S,15S)-epoxy-(5Z,8Z,11Z)-icosatrienoic acid
        "C(CCC/C=C\\C[C@@H]1/C(/O1)=C/C=C\\C/C=C\\CCCCC)(=O)O",  # (5Z,8R,9Z,11Z,14Z)-8,9-epoxyicosatetraenoic acid
    ]
    for s in test_smiles_list:
        decision, reason = is_cyclic_fatty_acid(s)
        print(f"SMILES: {s}")
        print(f"Result: {decision}")
        print(f"Reason: {reason}\n")