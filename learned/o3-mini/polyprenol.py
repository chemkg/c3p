"""
Classifies: CHEBI:26199 polyprenol
"""
#!/usr/bin/env python
"""
Classifies: polyprenol – any member of the class of prenols possessing the general formula
H-[CH2C(Me)=CHCH2]nOH in which n>=2 (i.e. more than one isoprene-like unit).

This improved classifier requires that:
  • The molecule is built exclusively from C, H and O.
  • The molecule is acyclic (a linear isoprenoid backbone) 
    (note: molecules with rings are not accepted even if they contain a long chain).
  • No carbonyl (C=O) groups are present.
  • It contains at least one terminal hydroxyl (-OH) group (on a primary carbon).
  • The longest carbon-only chain is sufficiently long (at least 10 C atoms).
  • The number of isoprene-like units (as judged via SMARTS patterns) is roughly what is expected
    from the length of that carbon chain. (A tolerance of ±1 unit is allowed.)
    
We define three SMARTS fragments for isoprene-like units:
    • alpha unit (terminal at the –OH side): [OH]CC([CH3])=C
    • internal unit: [CH2]C([CH3])=C[CH2]
    • omega unit (terminal at the other end): C([CH3])=C[CH2]O
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_polyprenol(smiles: str):
    """
    Determines if a molecule is a polyprenol based on its SMILES string.
    
    A polyprenol should:
       • Be composed solely of C, H and O.
       • Be acyclic.
       • Have no carbonyl (C=O) groups.
       • Possess at least one terminal hydroxyl (-OH) group on a primary carbon.
       • Possess a long carbon chain (at least 10 carbons in the longest chain).
       • Possess a number of isoprene-like repeating units that agrees with the chain length.
       
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if classified as a polyprenol, False otherwise.
        str: Explanation for the decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Check that molecule is solely built from C, H and O.
    allowed_atomic_nums = {1, 6, 8}
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in allowed_atomic_nums:
            return False, f"Molecule contains atom {atom.GetSymbol()} not allowed (only C, H, O are permitted)."
    
    # 2. Ensure molecule is acyclic.
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule contains rings; expected a linear isoprenoid chain."
    
    # 3. Exclude molecules with carbonyl groups (e.g. acids, esters)
    carbonyl_pat = Chem.MolFromSmarts("[CX3]=[OX1]")
    if mol.HasSubstructMatch(carbonyl_pat):
        return False, "Molecule contains carbonyl (C=O) groups which are not expected for a polyprenol."
    
    # 4. Check for at least one hydroxyl (-OH) group.
    hydroxyl_pat = Chem.MolFromSmarts("[OX2H]")
    if not mol.HasSubstructMatch(hydroxyl_pat):
        return False, "No hydroxyl (-OH) group found."

    # 5. Look for a terminal -OH on a primary carbon.
    terminal_oh_found = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 8:  # Oxygen
            # Must have at least one hydrogen -- using GetTotalNumHs.
            if atom.GetTotalNumHs() >= 1:
                # Check that oxygen is attached to exactly one heavy (non-hydrogen) atom.
                heavy_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() > 1]
                if len(heavy_neighbors) == 1:
                    neigh = heavy_neighbors[0]
                    if neigh.GetAtomicNum() == 6:
                        # Check that this carbon is primary: it must have only one carbon neighbor.
                        carbon_neighbors = [nbr for nbr in neigh.GetNeighbors() if nbr.GetAtomicNum() == 6]
                        if len(carbon_neighbors) == 1:
                            terminal_oh_found = True
                            break
    if not terminal_oh_found:
        return False, "No terminal hydroxyl (-OH) group on a primary carbon found."
    
    # 6. Check that there are a sufficient number of carbon atoms.
    carbon_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    n_carbons = len(carbon_atoms)
    if n_carbons < 10:
        return False, f"Too few carbons ({n_carbons}) for a polyprenol structure."
    
    # 7. Count isoprene-like repeating fragments via SMARTS.
    # Note: These SMARTS are heuristic.
    pattern_alpha = Chem.MolFromSmarts("[OH]CC([CH3])=C")
    pattern_internal = Chem.MolFromSmarts("[CH2]C([CH3])=C[CH2]")
    pattern_omega = Chem.MolFromSmarts("C([CH3])=C[CH2]O")
    
    matches_alpha = mol.GetSubstructMatches(pattern_alpha)
    matches_internal = mol.GetSubstructMatches(pattern_internal)
    matches_omega = mol.GetSubstructMatches(pattern_omega)
    total_units = len(matches_alpha) + len(matches_internal) + len(matches_omega)
    if total_units < 2:
        return False, f"Only {total_units} isoprene-like unit(s) detected (need at least 2)."
    
    # 8. Identify the longest carbon-only chain.
    # Build an undirected graph over carbon atoms only.
    carbon_indices = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    # Create an adjacency dictionary: key = carbon atom index, value = list of neighbor carbon indices.
    graph = {idx: [] for idx in carbon_indices}
    for idx in carbon_indices:
        atom = mol.GetAtomWithIdx(idx)
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6:
                graph[idx].append(nbr.GetIdx())
                
    # In a tree (acyclic) the longest chain (diameter) can be found by doing DFS from each terminal node.
    def dfs(start, visited):
        max_length = 1
        for nbr in graph[start]:
            if nbr not in visited:
                length = 1 + dfs(nbr, visited | {nbr})
                if length > max_length:
                    max_length = length
        return max_length

    # Terminal nodes in this carbon graph (degree 1).
    terminal_nodes = [idx for idx, nbrs in graph.items() if len(nbrs) == 1]
    longest_chain = 0
    for node in terminal_nodes:
        chain_length = dfs(node, {node})
        if chain_length > longest_chain:
            longest_chain = chain_length

    if longest_chain < 10:
        return False, f"Longest carbon chain has only {longest_chain} carbons; too short for a polyprenol."
    
    # 9. Estimate expected number of isoprene units.
    # In an ideal polyprenol of formula H-[CH2C(Me)=CHCH2]nOH, the chain has about n repeating 4-carbon units plus 1 extra carbon.
    expected_units = round((longest_chain - 1) / 4)
    # Allow a tolerance of ±1 unit.
    if abs(total_units - expected_units) > 1:
        return False, (f"Mismatch between repeating fragment count ({total_units}) and "
                       f"expected isoprene units (~{expected_units}) based on a {longest_chain}-carbon chain.")
    
    return True, (f"Polyprenol detected with {total_units} isoprene-like unit(s) "
                  f"and a {longest_chain}-carbon backbone (expected ~{expected_units} units), "
                  "with a terminal -OH group.")

# Example usage (you can remove or modify this block when using the function in another project):
if __name__ == "__main__":
    examples = {
        "(E,E,E)-geranylgeraniol": "CC(C)=CCC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C\\CO",
        "(2E,6E,10E)-omega-hydroxyfarnesol": "C\\C(CO)=C/CC\\C(C)=C\\CC\\C(C)=C\\CO",
        "Bionectin F": "OC(CCC=C(C)C)(CCCC(O)(CCCC(O)(CCCC(O)(CCCC(O)(CCCC(O)(CCCC(O)(CCC/C(=C/CC/C(=C/CO)/C)/C)C)C)C)C)C)C)C",
        "Dolichol-19": "OCC[C@H](CC/C=C(\\CC/C=C(\\CC/C=C(\\CC/C=C(\\CC/C=C(\\CC/C=C(\\CC/C=C(\\CC/C=C(\\CC/C=C(\\CC/C=C(\\CC/C=C(\\CC/C=C(\\CC\\C=C(\\CC/C=C(/CCC=C(C)C)\\C)/C)/C)/C)/C)/C)/C)/C)/C)/C)/C)/C)/C",
        "Glisoprenin E": "OC(C(O)CC/C(=C/CO)/C)(CC/C=C(/CC/C=C(/CCCC(O)(CCCC(O)(CCCC(O)(CCCC(O)(CCC(O)C(O)(C)C)C)C)C)C)\\C)\\C)C",
        "geraniol": "CC(C)=CCC\\C(C)=C\\CO",
        "Dolichol-18": "OCC[C@H](CC/C=C(\\CC/C=C(\\CC/C=C(\\CC/C=C(\\CC/C=C(\\CC/C=C(\\CC/C=C(\\CC/C=C(\\CC/C=C(\\CC/C=C(\\CC/C=C(\\CC/C=C(\\CC\\C=C(\\CC/C=C(/CCC=C(C)C)\\C)/C)/C)/C)/C)/C)/C)/C)/C)/C)/C)/C)/C)",
        "Glisoprenin A": "OC(CCC=C(C)C)(CCCC(O)(CCCC(O)(CCCC(O)(CCCC(=CCCC(=CCCC(=CCCC(=CCO)C)C)C)C)C)C)C)C",
        "(2-trans,6-trans)-farnesol": "CC(C)=CCC\\C(C)=C\\CC\\C(C)=C\\CO",
        "farnesol": "[H]C(CO)=C(C)CCC([H])=C(C)CCC=C(C)C",
        "SCH 66878": "OC(C(O)CC/C(=C/CC/C(=C/CC/C(=C/CO)/C)/C)/C)(CCCC(O)(CCCC(O)(CCCC(O)(CCCC(O)(CCCC(O)(CCCC(O)(CCCC(O)(CCCC(O)(CCCC(O)(CCC=C(C)C)C)C)C)C)C)C)C)C)C)C",
        "all-trans-octaprenol": "C(/C=C(/CC\\C=C(\\CC/C=C(\\C)/CCC=C(C)C)/C)\\C)C\\C(\\C)=C\\CC\\C(\\C)=C\\CC/C(/C)=C/CC\\C(=C\\CO)\\C",
        "ditrans,polycis-undecaprenol": "C(/C=C(\\CC/C=C(\\CC/C=C(\\CC/C=C(\\CC/C=C(\\CC/C=C(\\CC/C=C(\\CC/C=C(/CC/C=C(/CCC=C(C)C)\\C)\\C)/C)/C)/C)/C)/C)/C)/C)/C)O",
        "Gymnoprenol A10": "OC(C(O)CO)(CC/C=C(/CC/C=C(/CCCC(O)(CCCC(O)(CCCC(O)(CCCC(O)(CCCC(O)(CCCC(O)(CCC=C(C)C)C)C)C)C)C)C)\\C)\\C)C",
        "(2-cis,6-trans)-farnesol": "CC(C)=CCC\\C(C)=C\\CC\\C(C)=C/CO",
        "Glisoprenin D": "OC(CCC(O)C(O)(C)C)(CCCC(O)(CCCC(O)(CCCC(O)(CCC/C(=C/CC/C(=C/CC/C(=C/CC/C(=C/CO)/C)/C)/C)/C)C)C)C)C",
        "all-trans-undecaprenol": "CC(C)=CCC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C\\CO",
        "nerol": "C(=C\\CO)(\\CCC=C(C)C)/C",
        "(2-cis,6-cis)-farnesol": "CC(C)=CCC\\C(C)=C/CC\\C(C)=C/CO",
        "(2E,6E,10Z,14Z,18Z,22E)-3,7,11,15,19,23,27,31,35-Nonamethylhexatriaconta-2,6,10,14,18,22-hexaen-1-ol": "OC\\C=C(\\CC\\C=C(\\CC/C=C(\\CC/C=C(\\CC/C=C(\\CC\\C=C(\\CCCC(CCCC(CCCC(C)C)C)C)/C)/C)/C)/C)/C)/C",
        "geranylgeraniol": "CC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCO",
        "(2-trans,6-cis)-farnesol": "CC(C)=CCC\\C(C)=C/CC\\C(C)=C\\CO",
        "Glisoprenin F": "OC(CCC(O)C(O)(C)C)(CCCC(O)(CCCC(O)(CCCC(O)(CCCC(O)(CCC/C(=C/CC/C(=C/CC/C(=C/CO)/C)/C)/C)C)C)C)C)C",
        "all-trans-hexaprenol": "C(/C=C(/CCC=C(C)C)\\C)C\\C(\\C)=C\\CC\\C(\\C)=C\\CC/C(/C)=C/CC\\C(=C\\CO)\\C",
        "SCH 60057": "OC(CCC=C(C)C)(CCCC(O)(CCCC(O)(CCCC(O)(CCCC(O)(CCC/C(=C/CC/C(=C/CC/C(=C/CO)/C)/C)/C)C)C)C)C)C",
        "solanesol": "CC(C)=CCC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C\\CO"
    }
    
    for name, s in examples.items():
        result, explanation = is_polyprenol(s)
        print(f"{name}:\n  Result: {result}\n  Explanation: {explanation}\n")