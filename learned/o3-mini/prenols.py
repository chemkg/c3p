"""
Classifies: CHEBI:26244 prenols
"""
#!/usr/bin/env python3
"""
Classifies: Prenols
Definition: A prenol (or prenol derivative) is any acyclic alcohol possessing the general formula
  H-[CH2C(Me)=CHCH2]nOH with n>=1 – that is, its overall carbon count is 5 (for prenol) or a multiple
  of 5 (for longer polyprenyl chains), and (for molecules having >5 carbons) the isoprene-derived backbone
  contains at least one carbon–carbon double bond (and at least one methyl‐substituted double bond). In biosynthesis,
  the isoprene units are nearly always all‐trans. Terminal alcohols may be free or phosphorylated.
  
This implementation first verifies that the molecule is acyclic. It then extracts the carbon backbone
(as a tree) and computes the longest path in that subgraph so as to identify the chain ends. A valid prenol
must have the overall number of carbons equal to 5 or a multiple of 5 and one of the longest‐chain endpoints
must bear the terminal –OH (or –O–P) group. In addition, for molecules larger than prenol the presence of at least one
C=C double bond and at least one methyl‐substituted C=C (characteristic of the isoprene unit) is required, and any
double bond showing cis (Z) stereochemistry is rejected.
"""

from rdkit import Chem

def is_prenols(smiles: str):
    """
    Determines if a molecule belongs to the prenol class (or prenol derivative) based on its SMILES.
    Returns a Boolean value and a detailed reason.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Work with explicit hydrogens.
    mol = Chem.AddHs(mol)
    
    # 1. Must be acyclic.
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule contains rings which is atypical for prenols"
    
    # 2. Count total carbons.
    carbon_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    c_count = len(carbon_atoms)
    if c_count < 5:
        return False, f"Too few carbons ({c_count}) to be a prenol"
    if c_count > 5 and (c_count % 5 != 0):
        return False, f"Carbon count {c_count} is not a multiple of 5, which is atypical for prenols"
    
    # 3. Extract carbon backbone as a graph (only carbons and their C–C bonds).
    # Build a dictionary: key = carbon atom index, value = list of neighboring carbon atom indices.
    carbon_ids = set(atom.GetIdx() for atom in carbon_atoms)
    c_graph = { idx: [] for idx in carbon_ids }
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        if a1.GetIdx() in carbon_ids and a2.GetIdx() in carbon_ids:
            c_graph[a1.GetIdx()].append(a2.GetIdx())
            c_graph[a2.GetIdx()].append(a1.GetIdx())
    
    # If there are no bonds between carbons (should not happen for a chain), fail.
    if not any(c_graph.values()):
        return False, "No carbon–carbon bonds found in backbone"
    
    # 4. Compute the longest path in the carbon tree. For a tree, the longest path can be found via two DFS runs.
    def dfs(start, graph):
        # returns a dictionary of distances and a dictionary of predecessors from start.
        seen = {start: 0}
        prev = {start: None}
        stack = [start]
        while stack:
            node = stack.pop()
            for neigh in graph[node]:
                if neigh not in seen:
                    seen[neigh] = seen[node] + 1
                    prev[neigh] = node
                    stack.append(neigh)
        return seen, prev

    # Pick an arbitrary carbon from the graph.
    start = next(iter(carbon_ids))
    dist1, prev1 = dfs(start, c_graph)
    # Find farthest carbon from start.
    farthest = max(dist1, key=lambda k: dist1[k])
    dist2, prev2 = dfs(farthest, c_graph)
    # Now find the carbon farthest from 'farthest'.
    other_end = max(dist2, key=lambda k: dist2[k])
    # Reconstruct the longest path.
    path = []
    cur = other_end
    while cur is not None:
        path.append(cur)
        cur = prev2[cur]
    path = list(reversed(path))
    # The endpoints of the longest carbon chain:
    endpoints = {path[0], path[-1]}
    
    # 5. Check for terminal alcohol.
    # We search for free terminal –OH: SMARTS "[CH2][OX2H]" and phosphorylated terminal: "[CH2]O[P]"
    terminal_alc_smarts = Chem.MolFromSmarts("[CH2][OX2H]")
    terminal_alc_phos_smarts = Chem.MolFromSmarts("[CH2]O[P]")
    
    def is_carbon_terminal(carbon_idx):
        # Check that the carbon atom (with explicit H) has exactly one neighboring carbon.
        atom = mol.GetAtomWithIdx(carbon_idx)
        count_c = sum(1 for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 6)
        return (count_c == 1)
    
    has_terminal = False
    term_reason = ""
    # First check free –OH groups.
    matches_free = mol.GetSubstructMatches(terminal_alc_smarts)
    for match in matches_free:
        # match is a tuple (carbon_idx, oxygen_idx)
        c_idx = match[0]
        if c_idx in endpoints and is_carbon_terminal(c_idx):
            has_terminal = True
            break
    # Then check phosphorylated terminal alcohol.
    if not has_terminal:
        matches_phos = mol.GetSubstructMatches(terminal_alc_phos_smarts)
        for match in matches_phos:
            c_idx = match[0]
            if c_idx in endpoints and is_carbon_terminal(c_idx):
                has_terminal = True
                break
    if not has_terminal:
        return False, "No terminal alcohol (free or phosphorylated) found at a chain terminus"
    
    # 6. Check for at least one C=C double bond.
    dbl_bond_smarts = Chem.MolFromSmarts("[C]=[C]")
    if not mol.HasSubstructMatch(dbl_bond_smarts):
        return False, "No carbon–carbon double bond found, required for an isoprenoid structure"
    
    # 7. For molecules with more than 5 carbons, require at least one methyl‐substituted double bond.
    if c_count > 5:
        methyl_dbl_smarts = Chem.MolFromSmarts("[C]([CH3])=[C]")
        if not mol.HasSubstructMatch(methyl_dbl_smarts):
            return False, "No methyl‐substituted double bond ([C]([CH3])=[C]) detected, not typical of isoprene units"
    
    # 8. Check that all C=C bonds (between carbons) in the molecule are trans.
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
            a1 = bond.GetBeginAtom()
            a2 = bond.GetEndAtom()
            # Only consider bonds connecting two carbons.
            if a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 6:
                stereo = bond.GetStereo()
                # If stereo is defined and is Z (cis), then reject.
                if stereo == Chem.rdchem.BondStereo.STEREOZ:
                    return False, "Found a cis (Z) double bond; prenols are typically all‐trans"
    
    # 9. Compute expected number of isoprene units.
    if c_count == 5:
        expected_units = 1
    else:
        expected_units = (c_count // 5) - 1

    return True, f"Found terminal alcohol at a chain end; acyclic molecule with {c_count} carbons, consistent with {expected_units} isoprene repeating unit(s)."

# Testing examples when run as a script.
if __name__ == '__main__':
    test_smiles = [
        # True positives:
        "CC(C)=CCC\\C(C)=C\\CO",  # (E,E,E)-geranylgeraniol, 20 C => 3 repeating units.
        "C\\C(CO)=C/CC\\C(C)=C\\CC\\C(C)=C\\CO",  # (2E,6E,10E)-omega-hydroxyfarnesol, 15 C => 2 units.
        "CC(C)=CCC\\C(C)=C\\COP([O-])(=O)OP([O-])([O-])=O",  # hexaprenyl diphosphate(3-), 30 C.
        "CC(C)=CCC\\C(C)=C\\CO",  # geraniol, 10 C.
        "CC(C)=CCO",  # prenol, 5 C.
        "C\\C(CO)=C/CC\\C(C)=C\\C=O",  # (6E)-8-hydroxygeranial, 10 C.
        # False positives (should be rejected by some check – e.g. due to cis geometry, extra functional groups, etc.)
        "OC/C=C(\\CCC[C@@H](CCC[C@@H](CCCC(C)C)C)/C)/C",  # cis-Phytol (has a cis double bond)
        "CC(=O)C(O)C(=O)COP([O-])([O-])=O",  # 3-hydroxy-2,4-dioxopentyl phosphate(2-)
    ]
    for smi in test_smiles:
        result, reason = is_prenols(smi)
        print(f"SMILES: {smi}\n  Is prenol? {result}\n  Reason: {reason}\n")