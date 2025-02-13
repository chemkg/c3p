"""
Classifies: CHEBI:35915 sterol ester
"""
"""
Classifies: A sterol ester, defined as 
'A steroid ester obtained by formal condensation of the carboxy group of any carboxylic acid 
with the 3-hydroxy group of a sterol.' 

Detection strategy:
1. Look for at least one ester function (–O–C(=O)–).
2. Analyze the ring system to see if there is a fused ring system composed of exactly 4 rings,
   with one 5-membered ring and three 6-membered rings (i.e. the classical steroid skeleton).
3. Verify that one of the ester groups has its “O–” attached to an atom that is part of this 
   steroid nucleus.
"""

from rdkit import Chem

def is_sterol_ester(smiles: str):
    """
    Determines if a molecule is a sterol ester based on its SMILES string.
    The molecule should contain an ester group (O-C(=O)-) and a steroid nucleus 
    (four fused rings: three 6-membered and one 5-membered) with the ester oxygen 
    attached to the nucleus (the esterification of the 3-hydroxy group).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sterol ester, False otherwise
        str: Reason for classification
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # --- Step 1. Find ester groups ---
    # The SMARTS below matches an ester substructure O-C(=O)
    ester_smarts = "[O]-[C](=O)"
    ester_query = Chem.MolFromSmarts(ester_smarts)
    ester_matches = mol.GetSubstructMatches(ester_query)
    if not ester_matches:
        return False, "No ester group (O-C(=O)) found"
    
    # --- Step 2. Identify steroid nucleus via fused ring system ---
    ring_info = mol.GetRingInfo()
    # Get all rings as tuples (each ring is a tuple of atom indices)
    rings = list(ring_info.AtomRings())
    if not rings:
        return False, "No rings detected"
    
    # Build a simple graph among rings: two rings are "fused" if they share 2 or more atoms.
    num_rings = len(rings)
    adj = {i: set() for i in range(num_rings)}
    for i in range(num_rings):
        for j in range(i+1, num_rings):
            # Check if rings[i] and rings[j] share at least 2 atoms
            if len(set(rings[i]).intersection(rings[j])) >= 2:
                adj[i].add(j)
                adj[j].add(i)
    
    # Find connected components in the ring graph.
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
                for neigh in adj[node]:
                    if neigh not in comp:
                        stack.append(neigh)
            visited |= comp
            fused_components.append(comp)
    
    # Look for a fused ring component corresponding to the steroid nucleus:
    # Typically, the steroid nucleus has exactly 4 rings: one 5-membered and three 6-membered.
    steroid_atoms = set()
    steroid_found = False
    for comp in fused_components:
        if len(comp) != 4:
            continue  # not the 4-ring system we expect
        # Count ring sizes in this component.
        count5 = 0
        count6 = 0
        comp_atoms = set()
        for idx in comp:
            ring = rings[idx]
            comp_atoms.update(ring)
            if len(ring) == 5:
                count5 += 1
            elif len(ring) == 6:
                count6 += 1
        if count5 == 1 and count6 == 3:
            steroid_found = True
            steroid_atoms = comp_atoms
            break
    
    if not steroid_found:
        return False, "No steroid nucleus (4 fused rings with 1 five-membered and 3 six-membered) found"
    
    # --- Step 3. Check that one of the ester groups is attached to the steroid nucleus ---
    # In a sterol ester, the oxygen of the ester (that is not in the carbonyl) should be 
    # connected to an atom that is part of the steroid nucleus.
    for match in ester_matches:
        # The match tuple follows the order of atoms in the SMARTS pattern:
        # match[0] is the O, match[1] is the C of C(=O).
        o_idx = match[0]
        # Get the oxygen atom and its neighbors.
        o_atom = mol.GetAtomWithIdx(o_idx)
        neighbors = [nbr.GetIdx() for nbr in o_atom.GetNeighbors()]
        # In an ester O, one neighbor is the carbonyl carbon (the match) so we consider the other
        non_carbonyl_neighbors = [n for n in neighbors if n != match[1]]
        # Check if any non-carbonyl neighbor is part of the steroid nucleus.
        for n_idx in non_carbonyl_neighbors:
            if n_idx in steroid_atoms:
                return True, "Ester group is attached to a steroid nucleus (4 fused rings: 1 five- and 3 six-membered rings)"
    
    return False, "Ester group not attached to steroid nucleus via the expected 3-hydroxy position"

# For testing purposes you can try:
if __name__ == '__main__':
    # Example: cholesteryl arachidonate (one of the provided examples)
    test_smiles = "CCCCCC\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCC(=O)O[C@H]1CC[C@]2(C)[C@H]3CC[C@]4(C)[C@H](CC[C@H]4[C@@H]3CC=C2C1)[C@H](C)CCCC(C)C"
    result, reason = is_sterol_ester(test_smiles)
    print("Result:", result)
    print("Reason:", reason)