"""
Classifies: CHEBI:35915 sterol ester
"""
"""
Classifies: A sterol ester, defined as 
'A steroid ester obtained by formal condensation of the carboxy group of any carboxylic acid 
with the 3-hydroxy group of a sterol.'

Detection strategy:
1. Look for at least one ester function defined as an O connected to a carbonyl (O-C(=O)).
2. Detect a steroid nucleus by identifying one fused ring system containing exactly 4 rings:
   one 5-membered ring and three 6-membered rings.
3. Check that at least one ester group has its non-carbonyl oxygen directly attached to an atom
   within the steroid nucleus.
"""

from rdkit import Chem

def is_sterol_ester(smiles: str):
    """
    Determines if a molecule is a sterol ester based on its SMILES string.
    The molecule should contain an ester group (O-C(=O)-) and a steroid nucleus 
    (a fused system with exactly 4 rings, one 5-membered and three 6-membered)
    where the ester oxygen (the one not in the carbonyl) is attached directly to the nucleus.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a sterol ester, False otherwise
        str: Reason for classification
    """
    
    # Parse molecule from SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # --- Step 1. Find ester group(s) ---
    # This SMARTS matches an oxygen bonded to a carbon which is doubly bonded to oxygen.
    # Note: This will match an ester group O-C(=O) where the ester oxygen may have another neighbor.
    ester_smarts = "[O;D2]-[C](=O)"
    ester_query = Chem.MolFromSmarts(ester_smarts)
    ester_matches = mol.GetSubstructMatches(ester_query)
    if not ester_matches:
        return False, "No ester group (O-C(=O)) found"
    
    # --- Step 2. Identify steroid nucleus via fused ring system ---
    ring_info = mol.GetRingInfo()
    rings = list(ring_info.AtomRings())
    if not rings:
        return False, "No rings detected in molecule"
    
    # Build a simple graph between rings: two rings are fused if they share at least two atoms.
    num_rings = len(rings)
    adj = {i: set() for i in range(num_rings)}
    for i in range(num_rings):
        for j in range(i + 1, num_rings):
            if len(set(rings[i]).intersection(rings[j])) >= 2:
                adj[i].add(j)
                adj[j].add(i)
    
    # Find connected components (fused ring systems) from the ring graph.
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
    
    # Look for a fused ring component that corresponds to a steroid nucleus:
    # exactly 4 rings with one having 5 atoms and three having 6 atoms.
    steroid_atoms = set()
    steroid_found = False
    for comp in fused_components:
        if len(comp) != 4:
            continue  # Not the 4-ring system expected.
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
        return False, "No steroid nucleus found (expecting 4 fused rings: 1 five-membered and 3 six-membered)"
    
    # --- Step 3. Check ester attachment to the steroid nucleus ---
    # For a sterol ester, one of the ester oxygen atoms (the one not involved in the C=O)
    # must be directly attached to an atom of the steroid nucleus.
    for match in ester_matches:
        # In the match tuple, match[0] is the ester oxygen, match[1] is the carbonyl carbon.
        o_idx = match[0]
        o_atom = mol.GetAtomWithIdx(o_idx)
        # Get neighbors of the oxygen atom. One neighbor is the carbonyl carbon.
        neighbor_idxs = [nbr.GetIdx() for nbr in o_atom.GetNeighbors()]
        non_carbonyl_neighbors = [n for n in neighbor_idxs if n != match[1]]
        for n_idx in non_carbonyl_neighbors:
            if n_idx in steroid_atoms:
                return True, "Ester group is attached to a steroid nucleus (4 fused rings: 1 five- and 3 six-membered)"
    
    return False, "Ester group not attached to the steroid nucleus via the expected 3-hydroxy position"

# For testing purposes, you can try one of the provided examples:
if __name__ == '__main__':
    # Example: cholesteryl arachidonate
    test_smiles = "CCCCCC\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCC(=O)O[C@H]1CC[C@]2(C)[C@H]3CC[C@]4(C)[C@H](CC[C@H]4[C@@H]3CC=C2C1)[C@H](C)CCCC(C)C"
    result, reason = is_sterol_ester(test_smiles)
    print("Result:", result)
    print("Reason:", reason)