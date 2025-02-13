"""
Classifies: CHEBI:26208 polyunsaturated fatty acid
"""
"""
Classifies: Polyunsaturated fatty acid (PUFA)
Definition: Any fatty acid containing more than one double bond.
Acids in this group are reported to have cardioprotective effects;
and levels are lowered in chronic fatigue syndrome.

This implementation first checks for a free (terminal) carboxyl group.
We require that a carboxyl group is “free” – that is, the carboxyl carbon
has exactly one carbon neighbor (the beginning of the aliphatic chain).
Then we count carbon–carbon double bonds (ignoring the carbonyl C=O bond).
We also check that the molecule is acyclic and has a sufficiently long carbon chain.
"""

from rdkit import Chem

def is_polyunsaturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a polyunsaturated fatty acid based on its SMILES string.
    A polyunsaturated fatty acid (PUFA) must have a free carboxyl group (terminal acid),
    contain at least 2 carbon-carbon double bonds, be acyclic, and feature an aliphatic chain
    of sufficient length.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polyunsaturated fatty acid, False otherwise
        str: Explanation of the classification
    """
    # Parse SMILES to an RDKit molecule object.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define two SMARTS patterns for a free carboxyl group.
    # (Some acids under physiological conditions are deprotonated.)
    carboxyl_pattern1 = Chem.MolFromSmarts("[CX3](=O)[OX1H]")  # protonated acid
    carboxyl_pattern2 = Chem.MolFromSmarts("[CX3](=O)[O-]")    # deprotonated acid

    matches1 = mol.GetSubstructMatches(carboxyl_pattern1)
    matches2 = mol.GetSubstructMatches(carboxyl_pattern2)
    free_carboxyl_matches = matches1 + matches2

    if not free_carboxyl_matches:
        return False, "No free carboxyl group found"

    # We require that the carboxyl group be terminal (i.e. the carboxyl carbon has exactly one carbon neighbor).
    terminal_match = None
    for match in free_carboxyl_matches:
        # In our SMARTS the first atom is the carboxyl carbon.
        carboxyl_c = mol.GetAtomWithIdx(match[0])
        # Count the number of neighboring carbons (atomic num==6).
        carbon_neighbors = [nbr for nbr in carboxyl_c.GetNeighbors() if nbr.GetAtomicNum() == 6]
        if len(carbon_neighbors) == 1:
            terminal_match = match
            break
    if terminal_match is None:
        return False, "No terminal free carboxyl group found"
    
    # Identify the carbon neighbor that begins the fatty acid chain.
    carboxyl_c = mol.GetAtomWithIdx(terminal_match[0])
    chain_start = None
    for nbr in carboxyl_c.GetNeighbors():
        if nbr.GetAtomicNum() == 6:
            chain_start = nbr
            break
    if chain_start is None:
        return False, "Carboxyl carbon does not connect to any alkyl chain"
    
    # Count the number of carbon-carbon double bonds (ignoring the acid C=O bond).
    double_bond_count = 0
    for bond in mol.GetBonds():
        # Check only for C=C bonds.
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
            a1 = bond.GetBeginAtom()
            a2 = bond.GetEndAtom()
            if a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 6:
                double_bond_count += 1
    if double_bond_count < 2:
        return False, f"Only {double_bond_count} carbon-carbon double bond(s) found; need at least 2"
    
    # Reject molecules that contain rings as most free fatty acids are acyclic.
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule contains rings; not a typical free fatty acid"
    
    # Define a helper function: perform DFS on carbon atoms (ignoring non-carbon atoms)
    # to determine the longest contiguous carbon chain length.
    def dfs(atom, visited):
        max_length = 1
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6 and nbr.GetIdx() not in visited:
                visited.add(nbr.GetIdx())
                length = 1 + dfs(nbr, visited)
                visited.remove(nbr.GetIdx())
                if length > max_length:
                    max_length = length
        return max_length

    # Start the DFS from the chain_start. We count the chain length as the number of carbon atoms
    # in the longest contiguous (acyclic) carbon path. (We do not count the carboxyl carbon.)
    chain_length = dfs(chain_start, {chain_start.GetIdx()})
    
    # Require that the fatty acyl chain be long enough (threshold can be adjusted).
    # Here we require at least 4 carbons beyond the carboxyl carbon.
    if chain_length < 4:
        return False, f"Aliphatic chain too short (only {chain_length} carbons found); not a fatty acid"
    
    return True, f"Contains terminal carboxyl group, {double_bond_count} carbon-carbon double bonds, and an aliphatic chain of {chain_length} carbons, meets PUFA criteria"

# (Optional) simple tests; remove these lines if integrating into a larger package.
if __name__ == "__main__":
    test_smiles = [
        "CC(C)=CCCC(C)=CCCC(C)=CC(O)=O",  # farnesoic acid (should be True)
        "CCCCCCCC/C=C/CCCCCCCC(=O)OC[C@H](COP([O-])(=O)OC[C@H]([NH3+])C([O-])=O)OC(=O)CCC/C=C/C/C=C/C/C=C/C/C=C/CCCCC", # phospholipid false positive
        "OC(=O)CCCCCCC\\C=C\\C=C\\CCCCCC", # example with one double bond only (False)
    ]
    for smi in test_smiles:
        result, reason = is_polyunsaturated_fatty_acid(smi)
        print(f"SMILES: {smi}\nResult: {result}\nReason: {reason}\n")