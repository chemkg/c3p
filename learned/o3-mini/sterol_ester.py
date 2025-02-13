"""
Classifies: CHEBI:35915 sterol ester
"""
#!/usr/bin/env python
"""
Classifies: A sterol ester, defined as  
'A steroid ester obtained by formal condensation of the carboxy group of any carboxylic acid 
with the 3-hydroxy group of a sterol.'

Detection strategy:
1. Look for an ester group defined as an oxygen (with two bonds) attached to a carbonyl carbon.
2. Identify a steroid nucleus defined as a fused ring system composed of exactly four rings – one 5‐membered and three 6‐membered – 
   whose atoms are exclusively carbons and that share the expected connectivity (two rings must be terminal with only one linkage, etc.).
3. Confirm that for at least one ester group the “ester oxygen” (the oxygen not part of the carbonyl) is directly attached to an atom
   that lies within the steroid nucleus.
"""

from rdkit import Chem

def is_sterol_ester(smiles: str):
    """
    Determines if a molecule is a sterol ester based on its SMILES string.
    The molecule should contain an ester group (O-C(=O)-) and a steroid nucleus – defined here as a fused carbocyclic 
    system composed of exactly four rings (three six-membered and one five-membered) with proper connectivity—
    with at least one ester group attached (via its non-carbonyl oxygen) to that nucleus.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule is a sterol ester, False otherwise
        str: Reason for the classification
    """
    # --- Step 0. Parse molecule ---
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # --- Step 1. Locate ester groups ---
    # We are looking for the pattern: an oxygen (with 2 bonds) attached to a carbon which is double-bonded to oxygen.
    ester_smarts = "[O;D2]-[C](=O)"
    ester_query = Chem.MolFromSmarts(ester_smarts)
    ester_matches = mol.GetSubstructMatches(ester_query)
    if not ester_matches:
        return False, "No ester group (O-C(=O)) found"
    
    # --- Step 2. Identify steroid nucleus ---
    # Get all rings in the molecule
    ring_info = mol.GetRingInfo()
    rings = list(ring_info.AtomRings())
    if not rings:
        return False, "No rings detected in molecule"
    
    # Build ring adjacency: two rings are “fused” if they share at least 2 atoms.
    num_rings = len(rings)
    ring_adj = {i: set() for i in range(num_rings)}
    for i in range(num_rings):
        for j in range(i + 1, num_rings):
            if len(set(rings[i]).intersection(rings[j])) >= 2:
                ring_adj[i].add(j)
                ring_adj[j].add(i)
    
    # Find fused ring systems (connected components in ring_adj)
    visited = set()
    fused_components = []
    for i in range(num_rings):
        if i not in visited:
            comp = set()
            stack = [i]
            while stack:
                node = stack.pop()
                if node in comp:
                    continue
                comp.add(node)
                for neigh in ring_adj[node]:
                    if neigh not in comp:
                        stack.append(neigh)
            visited |= comp
            fused_components.append(comp)
    
    steroid_atoms = set()
    steroid_found = False
    for comp in fused_components:
        # We require exactly 4 fused rings.
        if len(comp) != 4:
            continue
        comp = list(comp)
        # Check ring sizes: we expect one 5-membered ring and three 6-membered rings.
        sizes = [len(rings[i]) for i in comp]
        if sorted(sizes) != [5, 6, 6, 6]:
            continue
        # Additionally, require that all atoms in these rings are carbons.
        comp_atom_set = set()
        for i in comp:
            comp_atom_set.update(rings[i])
        if any(mol.GetAtomWithIdx(idx).GetAtomicNum() != 6 for idx in comp_atom_set):
            continue
        # Build the subgraph of rings
        comp_adj = {}
        for i in comp:
            comp_adj[i] = set()
            for j in comp:
                if i != j and len(set(rings[i]).intersection(rings[j])) >= 2:
                    comp_adj[i].add(j)
        # Check connectivity: for a typical steroid nucleus we expect exactly two rings
        # with one connection (degree 1) and two rings with degree 2.
        degrees = {i: len(comp_adj[i]) for i in comp}
        deg_values = sorted(degrees.values())
        if deg_values != [1, 1, 2, 2]:
            continue
        # Additionally, require that the unique 5-membered ring (the A- or D-ring typically) is terminal (degree 1)
        five_ring_idxs = [i for i in comp if len(rings[i]) == 5]
        if not five_ring_idxs or any(degrees[i] != 1 for i in five_ring_idxs):
            continue
        
        # If we reach here, we consider this fused system our steroid nucleus.
        steroid_found = True
        steroid_atoms = comp_atom_set  # atoms that belong to the nucleus
        break

    if not steroid_found:
        return False, "No steroid nucleus found (expected 4 fused carbocyclic rings: 1 five-membered and 3 six-membered with proper connectivity)"
    
    # --- Step 3. Check that an ester group is attached to the steroid nucleus ---
    # For each ester match, the queried pattern gives us:
    #   match[0]: the non-carbonyl oxygen; match[1]: the carbonyl carbon.
    for match in ester_matches:
        o_idx, c_idx = match[0], match[1]
        o_atom = mol.GetAtomWithIdx(o_idx)
        # Look at all neighbors of the oxygen aside from the carbonyl carbon.
        neighbor_idxs = [nbr.GetIdx() for nbr in o_atom.GetNeighbors()]
        non_carbonyl_neighbors = [n for n in neighbor_idxs if n != c_idx]
        # If any such neighbor is part of the steroid nucleus then we count it as attached at the 3-hydroxy position.
        for n_idx in non_carbonyl_neighbors:
            if n_idx in steroid_atoms:
                return True, "Ester group is attached to a steroid nucleus via the expected 3-hydroxy position"
    return False, "No ester group is attached to the steroid nucleus via the expected position"

# For testing purposes (you may adjust or add examples):
if __name__ == '__main__':
    # Example: cholesteryl arachidonate (one known sterol ester)
    test_smiles = "CCCCCC\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCC(=O)O[C@H]1CC[C@]2(C)[C@H]3CC[C@]4(C)[C@H](CC[C@H]4[C@@H]3CC=C2C1)[C@H](C)CCCC(C)C"
    result, reason = is_sterol_ester(test_smiles)
    print("Result:", result)
    print("Reason:", reason)