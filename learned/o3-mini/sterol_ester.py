"""
Classifies: CHEBI:35915 sterol ester
"""
#!/usr/bin/env python
"""
Classifies: A sterol ester, defined as 
'A steroid ester obtained by formal condensation of the carboxy group of any carboxylic acid 
with the 3-hydroxy group of a sterol.'

Detection strategy:
1. Look for at least one ester function defined as an oxygen bonded to a carbon that is doubly bonded to oxygen.
2. Identify a steroid nucleus defined as a fused ring system composed of exactly four rings – one 5‐membered ring and three 6‐membered rings – 
   with the expected connectivity pattern (the two “terminal” rings have degree 1 and the 5‐membered ring must be terminal).
3. Confirm that for at least one ester group the “ester oxygen” (the oxygen not part of the carbonyl group) is directly attached to an atom 
   that lies within the steroid nucleus.
"""

from rdkit import Chem

def is_sterol_ester(smiles: str):
    """
    Determines if a molecule is a sterol ester based on its SMILES string.
    The molecule should contain an ester group (O-C(=O)-) and a steroid nucleus 
    (a fused system with exactly four rings, one 5-membered and three 6-membered, with the proper connectivity)
    and at least one ester group must be attached (through its non-carbonyl oxygen) to that nucleus.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a sterol ester, False otherwise.
        str: Reason for classification.
    """
    # --- Step 0. Parse molecule ---
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # --- Step 1. Locate ester groups ---
    # We look for pattern: an oxygen (with two bonds) attached to a carbonyl carbon.
    ester_smarts = "[O;D2]-[C](=O)"
    ester_query = Chem.MolFromSmarts(ester_smarts)
    ester_matches = mol.GetSubstructMatches(ester_query)
    if not ester_matches:
        return False, "No ester group (O-C(=O)) found"
    
    # --- Step 2. Identify steroid nucleus ---
    # First, get all rings in the molecule.
    ring_info = mol.GetRingInfo()
    rings = list(ring_info.AtomRings())
    if not rings:
        return False, "No rings detected in molecule"
    
    # Build a ring adjacency graph: two rings are connected if they share at least 2 atoms.
    num_rings = len(rings)
    ring_adj = {i: set() for i in range(num_rings)}
    for i in range(num_rings):
        for j in range(i + 1, num_rings):
            if len(set(rings[i]).intersection(rings[j])) >= 2:
                ring_adj[i].add(j)
                ring_adj[j].add(i)
    
    # Find connected components in the ring graph (each is a fused ring system).
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
    # Look into each fused component and require:
    # (a) It contains exactly 4 rings.
    # (b) The sizes of the rings are one 5-membered and three 6-membered.
    # (c) The pattern of shared fusions within the component is compatible with a steroid:
    #     - Build a graph (nodes = rings in comp, edges = fusion if share >=2 atoms).
    #     - In a typical steroid nucleus, the two terminal rings have degree 1,
    #       the two inner rings have degree 2, and the unique 5-membered ring (ring D) must have degree 1.
    for comp in fused_components:
        if len(comp) != 4:
            continue
        comp = list(comp)
        sizes = [len(rings[i]) for i in comp]
        if sorted(sizes) != [5,6,6,6]:
            continue

        # Build subgraph of the rings in this component.
        comp_adj = {}
        for i in comp:
            comp_adj[i] = set()
            for j in comp:
                if i != j and len(set(rings[i]).intersection(rings[j])) >= 2:
                    comp_adj[i].add(j)
        # Get degrees for rings in this component.
        degrees = {i: len(comp_adj[i]) for i in comp}
        # For a steroid nucleus, we expect exactly two rings with degree 1 and two rings with degree 2.
        degs = sorted(degrees.values())
        if degs != [1,1,2,2]:
            continue
        # Additionally, require that the unique 5-membered ring (if any) has degree 1.
        five_ring_idxs = [i for i in comp if len(rings[i]) == 5]
        if not five_ring_idxs:
            continue
        # In a proper steroid nucleus, the 5-membered ring should be terminal (degree 1).
        if any(degrees[i] != 1 for i in five_ring_idxs):
            continue
        # We now assume that this fused ring system is our steroid nucleus.
        steroid_found = True
        # Collect all atom indices that belong to any ring in this component.
        for i in comp:
            steroid_atoms.update(rings[i])
        break

    if not steroid_found:
        return False, "No steroid nucleus found (expected 4 fused rings: 1 five-membered and 3 six-membered with proper connectivity)"
    
    # --- Step 3. Check that an ester group is attached to the steroid nucleus ---
    # For each ester match the oxygen (match[0]) will be the non–carbonyl oxygen.
    for match in ester_matches:
        o_idx, c_idx = match[0], match[1]
        o_atom = mol.GetAtomWithIdx(o_idx)
        # Get the neighbors of the ester oxygen.
        neighbor_idxs = [nbr.GetIdx() for nbr in o_atom.GetNeighbors()]
        # Remove the carbonyl carbon.
        non_carbonyl = [n for n in neighbor_idxs if n != c_idx]
        # If any such neighbor belongs to the steroid nucleus, we consider it a valid attachment.
        for n_idx in non_carbonyl:
            if n_idx in steroid_atoms:
                return True, "Ester group is attached to a steroid nucleus via the expected 3-hydroxy position"
    return False, "Ester group not attached to the steroid nucleus via the expected position"

# For testing purposes (you may adjust or add examples):
if __name__ == '__main__':
    # Example: cholesteryl arachidonate
    test_smiles = "CCCCCC\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCC(=O)O[C@H]1CC[C@]2(C)[C@H]3CC[C@]4(C)[C@H](CC[C@H]4[C@@H]3CC=C2C1)[C@H](C)CCCC(C)C"
    result, reason = is_sterol_ester(test_smiles)
    print("Result:", result)
    print("Reason:", reason)