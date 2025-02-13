"""
Classifies: CHEBI:73155 trienoic fatty acid
"""
"""
Classifies: Trienoic fatty acid
Definition: Any polyunsaturated fatty acid (with only C, H, O) that contains exactly three carbon–carbon double bonds 
and has a free terminal carboxylic acid group. Additionally, the molecule should contain a sufficiently long carbon chain.
Heuristics:
  1. Molecule must contain only C, H, and O atoms.
  2. Must have a free (non‐ester) carboxylic acid group – the acid carbon must be terminal (attached to only one carbon) 
     and its immediate neighbor (alpha carbon) must not be over‐substituted (expected CH2).
  3. There must be exactly three pure C=C bonds (ignoring any carbonyl C=O bonds).
  4. The alkyl chain (starting at the alpha carbon) must be long enough (at least 8 carbon atoms) to be a fatty acid.
"""
from rdkit import Chem

def is_trienoic_fatty_acid(smiles: str):
    """
    Determines if a molecule is a trienoic fatty acid based on its SMILES string.
    A trienoic fatty acid (by our extended heuristics) must:
      • contain only C, H, and O atoms,
      • possess a free, terminal carboxylic acid group,
      • have exactly three pure carbon–carbon double bonds (C=C), and
      • include a sufficiently long carbon chain (at least 8 carbons in the chain beyond the acid group).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule meets our criteria, False otherwise.
        str: A reason message explaining the outcome.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Check that only atoms allowed (C, H, O) are present.
    allowed_atomic_nums = {1, 6, 8}  # H, C, O
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in allowed_atomic_nums:
            return False, "Contains atoms other than C, H, and O; not a fatty acid"
    
    # 2. Look for a free (non‐ester) carboxylic acid group.
    # We use two SMARTS to catch both protonated and deprotonated forms.
    acid_smarts_list = [
        "[CX3](=O)[OX2H1]",  # e.g. C(=O)O
        "[CX3](=O)[O-]"      # e.g. C(=O)[O-]
    ]
    acid_match = None
    for smarts in acid_smarts_list:
        acid_pat = Chem.MolFromSmarts(smarts)
        matches = mol.GetSubstructMatches(acid_pat)
        if matches:
            acid_match = matches[0]  # take the first found acid group
            break
    if acid_match is None:
        return False, "No carboxylic acid group found; not a fatty acid"

    # In our SMARTS the first atom is the acid carbon (the carbonyl carbon of the acid group).
    acidC_idx = acid_match[0]
    acidC = mol.GetAtomWithIdx(acidC_idx)
    # Check that the acid carbon is terminal – it should only be attached to one carbon.
    carbon_neighbors = [nbr for nbr in acidC.GetNeighbors() if nbr.GetAtomicNum() == 6]
    if len(carbon_neighbors) != 1:
        return False, "Carboxylic acid group is not terminal; may be esterified or not a free fatty acid"
    # Check the alpha carbon (neighboring carbon) for typical substitution.
    alphaC = carbon_neighbors[0]
    # For a canonical fatty acid, alphaC is usually CH2: expect 2 hydrogens.
    if alphaC.GetTotalNumHs() != 2:
        return False, "Alpha carbon adjacent to acid group is over/substituted; not a typical fatty acid backbone"
    
    # 3. Count the number of pure carbon–carbon double bonds.
    double_bond_count = 0
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
            a1 = bond.GetBeginAtom()
            a2 = bond.GetEndAtom()
            # Count only if both atoms are carbon (ignore C=O bonds)
            if a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 6:
                double_bond_count += 1
    if double_bond_count != 3:
        return False, f"Found {double_bond_count} carbon–carbon double bonds; need exactly 3"
    
    # 4. Check that the carbon chain is long enough.
    # We interpret a "fatty acid" as one having a long alkyl chain extending from the acid group.
    # We perform a DFS starting at the alpha carbon (only following carbon atoms) to find the longest chain.
    atom_indices = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    # Create an adjacency dictionary for carbon atoms.
    carbon_adj = {idx: [] for idx in atom_indices}
    for bond in mol.GetBonds():
        a_idx = bond.GetBeginAtomIdx()
        b_idx = bond.GetEndAtomIdx()
        if a_idx in carbon_adj and b_idx in carbon_adj:
            carbon_adj[a_idx].append(b_idx)
            carbon_adj[b_idx].append(a_idx)

    # DFS to find the maximum chain length from start_idx.
    def dfs(current, visited):
        max_length = 1  # count the current atom
        for neighbor in carbon_adj.get(current, []):
            if neighbor not in visited:
                length = 1 + dfs(neighbor, visited | {neighbor})
                if length > max_length:
                    max_length = length
        return max_length

    alpha_idx = alphaC.GetIdx()
    longest_chain = dfs(alpha_idx, {alpha_idx})
    # Heuristic: require at least 8 carbons in the chain (this threshold may be tuned).
    if longest_chain < 8:
        return False, f"Longest carbon chain starting from acid group has only {longest_chain} carbons; chain too short to be a fatty acid"
    
    # Passed all tests.
    return True, "Molecule is a trienoic fatty acid (free terminal carboxylic acid group, exactly three C=C bonds, and a sufficiently long alkyl chain)"

# Example testing; you can remove or comment these out in production.
if __name__ == "__main__":
    test_smiles = [
        "OC(=O)CCCCCCC/C=C\\CC(=O)/C=C(\\O)/C=C\\CCO",  # True positive
        "CCCC\\C=C/C\\C=C/CC(O)C(O)C\\C=C/CCCC([O-])=O",  # False positive example
        "CC\\C=C/C\\C=C/C\\C=C/CCCCCC[C@@H](OO)C(O)=O"     # False negative expected due to over-substituted alpha carbon
    ]
    for smi in test_smiles:
        result, reason = is_trienoic_fatty_acid(smi)
        print(f"SMILES: {smi}\n  Result: {result}\n  Reason: {reason}\n")