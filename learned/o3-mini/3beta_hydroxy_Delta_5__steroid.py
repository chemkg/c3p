"""
Classifies: CHEBI:1722 3beta-hydroxy-Delta(5)-steroid
"""
#!/usr/bin/env python
"""
Classifies: 3β-hydroxy-Δ(5)-steroid

A compound is considered to be a 3β-hydroxy-Δ(5)-steroid if it contains:
  1. A fused steroid nucleus – heuristically defined as a set of four fused rings (which may be a subset of a larger ring system)
     composed of exactly three six-membered rings and one five-membered ring.
  2. A hydroxyl group (–OH) attached to one of the carbons of that core fused ring set (approximating the 3β-OH group).
  3. At least one double bond between two carbons in the core fused system (approximating the Δ(5) double bond).

Note: The heuristic approach may miss some edge cases.
"""

from rdkit import Chem
from rdkit.Chem import AllChem
import itertools

def is_3beta_hydroxy_Delta_5__steroid(smiles: str):
    """
    Determines if a molecule is a 3β-hydroxy-Δ(5)-steroid based on heuristic structural features.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule fits the desired steroid class, False otherwise.
        str: Explanation for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return False, "Invalid SMILES string"
    
    # Get ring information as a list of tuples (each tuple: atom indices in that ring)
    ring_atoms = mol.GetRingInfo().AtomRings()
    if not ring_atoms:
        return False, "No rings found in molecule"
    
    n_rings = len(ring_atoms)
    # Build connectivity graph for rings:
    # Two rings are considered fused if they share >= 2 atoms.
    ring_graph = {i: set() for i in range(n_rings)}
    for i in range(n_rings):
        for j in range(i+1, n_rings):
            shared = set(ring_atoms[i]).intersection(ring_atoms[j])
            if len(shared) >= 2:
                ring_graph[i].add(j)
                ring_graph[j].add(i)
                
    # Find connected components among rings (each component is a candidate fused ring system)
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
        
    # Select the largest fused ring system (by union of atom indices)
    largest_comp = None
    largest_comp_atoms = set()
    comp_ring_list = []  # list of ring indices for each component
    for comp in components:
        atoms_in_comp = set()
        for rid in comp:
            atoms_in_comp.update(ring_atoms[rid])
        if len(atoms_in_comp) > len(largest_comp_atoms):
            largest_comp_atoms = atoms_in_comp
            largest_comp = comp

    # In the largest fused system, we have several rings. Instead of requiring exactly 4 rings,
    # we now try to find ANY subset of 4 rings (from the rings in the largest system) that meet the criteria:
    # exactly 3 rings of size 6 and 1 ring of size 5.
    candidate_rings = []
    for rid in largest_comp:
        candidate_rings.append(ring_atoms[rid])
    if len(candidate_rings) < 4:
        return False, "Largest fused ring system contains fewer than 4 rings"

    valid_subset = None
    for subset in itertools.combinations(candidate_rings, 4):
        six_membered = sum(1 for ring in subset if len(ring) == 6)
        five_membered = sum(1 for ring in subset if len(ring) == 5)
        if six_membered == 3 and five_membered == 1:
            valid_subset = subset
            break
    if not valid_subset:
        return False, ("No subset of 4 fused rings (3 six-membered and 1 five-membered) found in largest fused system; "
                       "found rings with sizes: " + ", ".join(str(len(r)) for r in candidate_rings))
    
    # Use the union of atoms in the valid subset as the steroid nucleus core.
    steroid_core_atoms = set()
    for ring in valid_subset:
        steroid_core_atoms.update(ring)
        
    # Look for a hydroxyl group (-OH) attached to a carbon in the steroid core.
    valid_hydroxyl_found = False
    for idx in steroid_core_atoms:
        atom = mol.GetAtomWithIdx(idx)
        if atom.GetAtomicNum() != 6:  # only check carbons
            continue
        # Count how many neighbors are also in the steroid core
        ring_neighbors = sum(1 for nbr in atom.GetNeighbors() if nbr.GetIdx() in steroid_core_atoms)
        # Check for an oxygen (that is bound to at least one hydrogen) attached to the carbon.
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 8:
                if nbr.GetTotalNumHs() >= 1:
                    # If the carbon is well embedded (at least 2 neighbors in core), then mark the hydroxyl as valid.
                    if ring_neighbors >= 2:
                        valid_hydroxyl_found = True
                        break
        if valid_hydroxyl_found:
            break
    if not valid_hydroxyl_found:
        return False, "No valid 3β-hydroxyl group found on a carbon of the steroid core fused system"
    
    # Look for a double bond (bond order == DOUBLE) within the steroid core.
    fused_double_bond_found = False
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.BondType.DOUBLE:
            a1 = bond.GetBeginAtomIdx()
            a2 = bond.GetEndAtomIdx()
            if a1 in steroid_core_atoms and a2 in steroid_core_atoms:
                fused_double_bond_found = True
                break
    if not fused_double_bond_found:
        return False, "No double bond found within the steroid core fused system (Δ(5) double bond expected)"
    
    return True, ("Molecule contains a steroid nucleus (3 six-membered and 1 five-membered rings in a fused core), "
                  "a 3β-hydroxyl group, and a Δ(5) double bond.")

# Example usage (if run as a script):
if __name__ == "__main__":
    # Example: Cholesterol is a typical 3β-hydroxy-Δ(5)-steroid.
    cholesterol_smiles = "C1[C@@]2([C@]3(CC[C@]4([C@]([C@@]3(CC=C2C[C@H](C1)O)[H])" \
                          "(CC[C@@]4([C@H](C)CCCC(C)C)[H])[H])C)[H])C"
    result, reason = is_3beta_hydroxy_Delta_5__steroid(cholesterol_smiles)
    print("Cholesterol classification:", result)
    print("Reason:", reason)