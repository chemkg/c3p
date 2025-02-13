"""
Classifies: CHEBI:1722 3beta-hydroxy-Delta(5)-steroid
"""
#!/usr/bin/env python
"""
Classifies: 3β-hydroxy-Δ(5)-steroid

A compound is considered to be a 3β-hydroxy-Δ(5)-steroid if it contains:
  1. A fused steroid nucleus – heuristically defined as a set of four fused rings (which may be a subset of a larger ring system)
     composed of exactly three six-membered rings and one five-membered ring. In addition, the union of the four rings
     should carry roughly 17 carbon atoms (allowing 16–18).
  2. A hydroxyl group (-OH) attached to one of the carbons in that nucleus (approximating the 3β-OH).
  3. At least one carbon–carbon double bond connecting two atoms of that nucleus (approximating the Δ(5) double bond).

Note: This heuristic is imperfect and may miss some edge cases.
"""

from rdkit import Chem
import itertools

def is_3beta_hydroxy_Delta_5__steroid(smiles: str):
    """
    Determines if a molecule is a 3β-hydroxy-Δ(5)-steroid based on a heuristic of its core fused ring system,
    the presence of a hydroxyl group on that core, and a double bond inside that core.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule is classified as a 3β-hydroxy-Δ(5)-steroid, False otherwise.
        str: Explanation for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Get ring information: a list of tuples (each tuple is the atom indices in that ring)
    ring_atoms = mol.GetRingInfo().AtomRings()
    if not ring_atoms:
        return False, "No rings found in molecule"
    
    n_rings = len(ring_atoms)
    # Build a graph of rings; two rings are considered fused if they share ≥ 2 atoms.
    ring_graph = {i: set() for i in range(n_rings)}
    for i in range(n_rings):
        for j in range(i+1, n_rings):
            if len(set(ring_atoms[i]).intersection(ring_atoms[j])) >= 2:
                ring_graph[i].add(j)
                ring_graph[j].add(i)
    
    # Find connected components among rings (each component is a candidate fused set of rings)
    visited = set()
    components = []
    for i in range(n_rings):
        if i not in visited:
            stack = [i]
            comp = set()
            while stack:
                cur = stack.pop()
                if cur not in visited:
                    visited.add(cur)
                    comp.add(cur)
                    stack.extend(ring_graph[cur] - visited)
            components.append(comp)
    if not components:
        return False, "No fused ring system found"

    # We now try to find within any connected component a subset of EXACTLY four rings 
    # that meets our criteria: exactly 3 rings of size 6 and 1 ring of size 5,
    # and the union of these ring atoms (the steroid nucleus) should contain roughly 17 carbon atoms.
    valid_subset = None
    valid_subset_atoms = None
    for comp in components:
        candidate_ring_sets = [ring_atoms[rid] for rid in comp]
        if len(candidate_ring_sets) < 4:
            continue
        # Try every combination of 4 rings from the component.
        for subset in itertools.combinations(candidate_ring_sets, 4):
            six_membered = sum(1 for ring in subset if len(ring) == 6)
            five_membered = sum(1 for ring in subset if len(ring) == 5)
            if six_membered == 3 and five_membered == 1:
                # Get the union of atoms in this candidate steroid nucleus.
                nucleus_atoms = set()
                for ring in subset:
                    nucleus_atoms.update(ring)
                # Count carbon atoms in candidate nucleus.
                c_count = 0
                for idx in nucleus_atoms:
                    atom = mol.GetAtomWithIdx(idx)
                    if atom.GetAtomicNum() == 6:
                        c_count += 1
                if c_count < 16 or c_count > 18:
                    # This subset does not have the expected carbon count for a steroid nucleus.
                    continue
                # If we reach here, we have a candidate nucleus.
                valid_subset = subset
                valid_subset_atoms = nucleus_atoms
                break
        if valid_subset is not None:
            break

    if valid_subset is None:
        # Provide detailed info based on available ring sizes in the largest fused system.
        all_sizes = []
        for comp in components:
            for rid in comp:
                all_sizes.append(len(ring_atoms[rid]))
        return False, ("No valid fused steroid nucleus found with 4 rings (3 six-membered and 1 five-membered) and ~17 carbons; "
                       "ring sizes available: " + ", ".join(str(sz) for sz in all_sizes))
    
    # Check for a hydroxyl group (–OH) attached to a carbon in the steroid nucleus.
    hydroxyl_found = False
    for idx in valid_subset_atoms:
        atom = mol.GetAtomWithIdx(idx)
        if atom.GetAtomicNum() != 6:  # only consider carbons in the nucleus
            continue
        # Look for an oxygen neighbor that is bound to at least one hydrogen (indicative of an OH)
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 8 and nbr.GetTotalNumHs() >= 1:
                # Assume that if this carbon is well-embedded (has ≥2 neighbors within nucleus), 
                # then the OH is part of the steroid core.
                core_neighbors = sum(1 for n in atom.GetNeighbors() if n.GetIdx() in valid_subset_atoms)
                if core_neighbors >= 2:
                    hydroxyl_found = True
                    break
        if hydroxyl_found:
            break
    if not hydroxyl_found:
        return False, "No valid 3β-hydroxyl group found attached to the steroid nucleus core"
    
    # Look for a double bond within the steroid nucleus.
    double_bond_found = False
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.BondType.DOUBLE:
            a1 = bond.GetBeginAtomIdx()
            a2 = bond.GetEndAtomIdx()
            if a1 in valid_subset_atoms and a2 in valid_subset_atoms:
                double_bond_found = True
                break
    if not double_bond_found:
        return False, "No double bond found in the steroid nucleus (Δ(5) double bond expected)"
    
    return True, ("Molecule contains a steroid nucleus (3 six-membered and 1 five-membered fused rings with ~17 carbons), "
                  "a 3β-hydroxyl group on the nucleus, and at least one double bond in the nucleus (Δ(5) bond).")

# Example usage (if run as a script):
if __name__ == "__main__":
    # Example: Cholesterol is a typical 3β-hydroxy-Δ(5)-steroid.
    cholesterol_smiles = "C1[C@@]2([C@]3(CC[C@]4([C@]([C@@]3(CC=C2C[C@H](C1)O)[H])(CC[C@@]4([C@H](C)CCCC(C)C)[H])[H])C)[H])C"
    result, reason = is_3beta_hydroxy_Delta_5__steroid(cholesterol_smiles)
    print("Cholesterol classification:", result)
    print("Reason:", reason)