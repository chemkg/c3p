"""
Classifies: CHEBI:1722 3beta-hydroxy-Delta(5)-steroid
"""
#!/usr/bin/env python
"""
Classifies: 3β-hydroxy-Δ(5)-steroid

A compound is considered to be a 3β-hydroxy-Δ(5)-steroid if it contains:
  1. A fused steroid nucleus – heuristically defined as a single (largest) fused ring system
     composed of exactly four rings, with three six-membered rings and one five-membered ring.
  2. A hydroxyl group (–OH) attached to a carbon in the fused system that is connected on at least 
     two sides to other atoms in the fused system (approximating the 3β-OH group).
  3. At least one double bond between two carbons that are both part of the fused system (approximating 
     the Δ(5) double bond).

Note: This heuristic approach may miss some edge cases or include some false positives.
"""

from rdkit import Chem
from rdkit.Chem import AllChem

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
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get ring information (list of tuples of atom indices in each ring)
    ring_atoms = mol.GetRingInfo().AtomRings()
    if not ring_atoms:
        return False, "No rings found"
    
    # Build connectivity graph among rings: two rings are connected in a fused system if they share at least 2 atoms
    n_rings = len(ring_atoms)
    ring_graph = {i: set() for i in range(n_rings)}
    for i in range(n_rings):
        for j in range(i+1, n_rings):
            shared = set(ring_atoms[i]).intersection(ring_atoms[j])
            if len(shared) >= 2:
                ring_graph[i].add(j)
                ring_graph[j].add(i)
    
    # Find connected components of rings (each component represents a candidate fused ring system)
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
    
    # Choose the largest fused ring system (by union of atom indices)
    largest_comp_atoms = set()
    largest_comp = None
    for comp in components:
        atoms_in_comp = set()
        for rid in comp:
            atoms_in_comp.update(ring_atoms[rid])
        if len(atoms_in_comp) > len(largest_comp_atoms):
            largest_comp_atoms = atoms_in_comp
            largest_comp = comp

    # Now count how many rings (from the candidate fused system) are entirely contained in the large fused system.
    # We require exactly 4 rings in the fused system with the required size:
    comp_rings = []
    for rid in largest_comp:
        ring = ring_atoms[rid]
        # Only consider rings all of whose atoms are in the fused union; they should be.
        if set(ring).issubset(largest_comp_atoms):
            comp_rings.append(ring)

    if len(comp_rings) != 4:
        return False, ("Fused ring system does not have exactly 4 rings (found {} rings); "
                       "steroid nucleus expected to be composed of 4 fused rings".format(len(comp_rings)))
    
    # Count six-membered and five-membered rings
    six_membered = sum(1 for ring in comp_rings if len(ring)==6)
    five_membered = sum(1 for ring in comp_rings if len(ring)==5)
    if six_membered != 3 or five_membered != 1:
        return False, ("Fused ring system does not match steroid nucleus ring sizes "
                       "(expected 3 six-membered and 1 five-membered, found {} six-membered and {} five-membered)"
                       .format(six_membered, five_membered))
    
    # Look for a hydroxyl group on a carbon that is part of the fused system.
    # For each carbon in the fused system, check if it has an –OH attached.
    valid_hydroxyl_found = False
    for idx in largest_comp_atoms:
        atom = mol.GetAtomWithIdx(idx)
        if atom.GetAtomicNum() != 6:  # only consider carbons
            continue
        # Count neighbors that are in the fused ring system
        ring_neighbors = sum(1 for nbr in atom.GetNeighbors() if nbr.GetIdx() in largest_comp_atoms)
        # Look for an –OH attached to this carbon
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 8:
                # Ensure the oxygen is in an OH group: it should have at least one attached hydrogen.
                # (Using GetTotalNumHs as a proxy.)
                if nbr.GetTotalNumHs() >= 1:
                    # We assume that if the carbon is part of the fused system and is connected to another fused atom (at least two),
                    # then an –OH here is likely the 3β-hydroxyl.
                    if ring_neighbors >= 2:
                        valid_hydroxyl_found = True
                        break
        if valid_hydroxyl_found:
            break
    if not valid_hydroxyl_found:
        return False, "No valid 3β-hydroxyl group found on a carbon of the fused steroid nucleus"
    
    # Look for a double bond (i.e. bond order = 2) inside the fused ring system.
    fused_double_bond_found = False
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.BondType.DOUBLE:
            a1 = bond.GetBeginAtomIdx()
            a2 = bond.GetEndAtomIdx()
            if a1 in largest_comp_atoms and a2 in largest_comp_atoms:
                fused_double_bond_found = True
                break
    if not fused_double_bond_found:
        return False, "No double bond found within the fused ring system (Δ(5) double bond expected)"
    
    return True, "Molecule contains a steroid nucleus (3 six-membered and 1 five-membered fused rings), a 3β-hydroxyl group, and a Δ(5) double bond"


# Example usage (if run as a script):
if __name__ == "__main__":
    # Cholesterol is a typical 3β-hydroxy-Δ(5)-steroid.
    cholesterol_smiles = "C1[C@@]2([C@]3(CC[C@]4([C@]([C@@]3(CC=C2C[C@H](C1)O)[H])" \
                          "(CC[C@@]4([C@H](C)CCCC(C)C)[H])[H])C)[H])C"
    result, reason = is_3beta_hydroxy_Delta_5__steroid(cholesterol_smiles)
    print("Cholesterol classification:", result)
    print("Reason:", reason)