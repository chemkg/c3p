"""
Classifies: CHEBI:39418 straight-chain saturated fatty acid
"""
#!/usr/bin/env python
"""
Classifies: CHEBI: Straight-chain saturated fatty acid (defined as any saturated fatty acid lacking a side‐chain)

A straight-chain fatty acid must contain a single carboxylic acid group and an unbranched carbon skeleton
(i.e. when considering only carbon–carbon single bonds, the carbon atoms form a single simple “path”
with two endpoints). Hydroxyl or oxo modifications attached to the chain (which are not carbons) are allowed.
Molecules with extra carbon substituents (side-chains) or any carbon–carbon unsaturation will be rejected.
"""
from rdkit import Chem

def is_straight_chain_saturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a straight-chain saturated fatty acid (no carbon side-chains)
    based solely on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a straight-chain saturated fatty acid, False otherwise.
        str: Explanation of the result.
    """
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for carboxylic acid substructure.
    # First try to match the protonated acid (-C(=O)O) pattern.
    acid_smarts = Chem.MolFromSmarts("C(=O)[O;H]")
    acid_matches = mol.GetSubstructMatches(acid_smarts)
    if not acid_matches:
        # If no match, try deprotonated acid (-C(=O)[O-])
        acid_smarts2 = Chem.MolFromSmarts("C(=O)[O-]")
        acid_matches = mol.GetSubstructMatches(acid_smarts2)
        if not acid_matches:
            return False, "No carboxylic acid group found"

    # Assume the first match; the acid carbon is the first atom in the matched tuple.
    acid_C_idx = acid_matches[0][0]

    # Build a graph only from carbon atoms connected by single (C-C) bonds.
    # In a true straight-chain fatty acid, the carboxyl carbon and the alkyl chain
    # will form a single linear (path) component.
    carbon_indices = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    # Create an adjacency dictionary for all carbon atoms based on single bonds.
    adjacency = {idx: set() for idx in carbon_indices}
    for bond in mol.GetBonds():
        # Consider only single bonds for the carbon skeleton
        if bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
            a1 = bond.GetBeginAtomIdx()
            a2 = bond.GetEndAtomIdx()
            if a1 in adjacency and a2 in adjacency:
                adjacency[a1].add(a2)
                adjacency[a2].add(a1)

    if acid_C_idx not in adjacency:
        return False, "Carboxyl carbon not recognized as part of the carbon skeleton"

    # Find the connected component in the carbon-only graph that contains the acid carbon.
    component = set()
    stack = [acid_C_idx]
    while stack:
        current = stack.pop()
        if current not in component:
            component.add(current)
            stack.extend(adjacency[current] - component)

    # For a linear (unbranched) chain (i.e. a simple path), the graph should have exactly 2 endpoints 
    # (vertices with one neighbor) and all other vertices (if any) should have exactly 2 neighbors.
    n = len(component)
    if n < 2:
        return False, "Fatty acid chain too short"
    # Get degree for each vertex (only counting connections within the component)
    degree_counts = [len(adjacency[idx] & component) for idx in component]
    end_points = sum(1 for d in degree_counts if d == 1)
    middle_points = sum(1 for d in degree_counts if d == 2)
    # Any vertex with degree not equal to 1 or 2 signals branching or an abnormal chain.
    others = [d for d in degree_counts if d not in (1, 2)]
    
    if n == 2:
        # For two carbons the only possibility is that both are endpoints.
        if end_points != 2:
            return False, "Unexpected bonding in a minimal chain"
    else:
        if not (end_points == 2 and middle_points == (n - 2) and len(others) == 0):
            return False, "Branching or non-linear carbon chain detected"
            
    # At this point, the carbon skeleton (including the acid carbon) is linear.
    # Note that any non-carbon substituents (e.g. hydroxy, oxo, deuterium) do not affect our test.
    # Also, since we only considered C-C single bonds, any carbon–carbon unsaturation (like C=C)
    # would have been excluded from the path.
    return True, "Molecule is a straight-chain saturated fatty acid lacking side-chains"
    
# Uncomment the lines below to test examples:
# test_smiles = [
#     "CCCCCCCCCCCCCCCCCCCC(O)=O",  # icosanoic acid
#     "CCCCCCCCCCCCCCCCCCCCCCCCCCCCC(O)=O",  # hexacosanoic acid
#     "C[C@@H](O)CCCCCCCCCCCCCCCCCCC(O)=O",  # (20R)-20-hydroxyhenicosanoic acid
#     "OC(C)CCCCCCCCCCCCCCCCCCC(=O)O",  # 20-hydroxyhenicosanoic acid
#     "CCCC(O)=O",  # butyric acid
#     "C(CCCCCC)CCC(=O)O",  # decanoic acid with evident branching -> should fail
# ]
# for s in test_smiles:
#     res, reason = is_straight_chain_saturated_fatty_acid(s)
#     print(f"SMILES: {s}\nResult: {res}\nReason: {reason}\n")