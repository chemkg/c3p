"""
Classifies: CHEBI:35343 17beta-hydroxy steroid
"""
"""
Classifies: CHEBI:17beta-hydroxy steroid
Definition: A 17-hydroxy steroid in which the hydroxy group at position 17 has a beta-configuration.
Heuristic update:
  1. Find all carbocyclic rings (rings composed solely of carbon) of size 5 or 6.
  2. Build a graph among these rings where an edge connects two rings sharing at least 2 atoms.
  3. For each connected component, iterate over all subsets of 4 rings. Accept a subset only if it
     contains exactly one 5-membered ring and three 6-membered rings (the classical steroid nucleus).
  4. Within the selected steroid core (the union of atoms in that subset), check for a candidate 17beta-hydroxy group:
     a) A carbon atom (inside the core) that is chiral (has defined stereochemistry) and
     b) is a member of the unique 5-membered ring (the “D ring” candidate)
     c) and has at least one oxygen neighbor connected by a single bond.
     
If such a candidate is found, return True with reasoning; otherwise, return False.
"""

from rdkit import Chem
from itertools import combinations

def is_17beta_hydroxy_steroid(smiles: str):
    """
    Determines whether a molecule is a 17beta-hydroxy steroid.
    
    The algorithm finds four fused carbocyclic rings (three of 6 atoms and one of 5 atoms)
    and then in the unique 5-membered ring looks for a chiral carbon carrying an -OH group.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a 17beta-hydroxy steroid, False otherwise.
        str: Explanation for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Ensure chiral information is computed.
    Chem.AssignStereochemistry(mol, cleanIt=True)
    
    # Get all rings present in the molecule
    ring_info = mol.GetRingInfo().AtomRings()
    # Only keep carbocyclic rings (only carbon atoms) and that are 5- or 6-membered.
    carbocyclic_rings = []
    for ring in ring_info:
        if len(ring) in (5, 6):
            if all(mol.GetAtomWithIdx(idx).GetAtomicNum() == 6 for idx in ring):
                carbocyclic_rings.append(set(ring))
                
    if not carbocyclic_rings:
        return False, "No 5- or 6-membered carbocyclic rings found; unlikely to be a steroid."
    
    # Build graph: each ring is a node; connect two rings if they share at least 2 atoms.
    n = len(carbocyclic_rings)
    ring_graph = {i: set() for i in range(n)}
    for i in range(n):
        for j in range(i+1, n):
            if len(carbocyclic_rings[i].intersection(carbocyclic_rings[j])) >= 2:
                ring_graph[i].add(j)
                ring_graph[j].add(i)
                
    # Find connected components among these rings (by their indices).
    visited = set()
    components = []
    for i in range(n):
        if i not in visited:
            comp = set()
            stack = [i]
            while stack:
                node = stack.pop()
                if node not in comp:
                    comp.add(node)
                    stack.extend(ring_graph[node] - comp)
            visited |= comp
            components.append(comp)
    
    # For each connected set of rings, try to find a subset of exactly 4 rings (3 six-membered and 1 five-membered).
    steroid_core_atoms = None
    five_ring_indices = None
    subset_found = False
    for comp in components:
        if len(comp) < 4:
            continue
        comp_rings = [ (i, carbocyclic_rings[i]) for i in comp ]
        # Iterate over all combinations of 4 rings in this connected component.
        for subset in combinations(comp_rings, 4):
            # Count 5- and 6-membered rings in this subset.
            count_5 = sum(1 for (_, ring) in subset if len(ring) == 5)
            count_6 = sum(1 for (_, ring) in subset if len(ring) == 6)
            if count_5 == 1 and count_6 == 3:
                # We accept this combination as the steroid core candidate.
                core_atoms = set()
                candidate_five_ring = None
                for (idx, ring) in subset:
                    core_atoms = core_atoms.union(ring)
                    if len(ring) == 5:
                        candidate_five_ring = ring
                steroid_core_atoms = core_atoms
                five_ring_indices = candidate_five_ring
                subset_found = True
                break
        if subset_found:
            break

    if steroid_core_atoms is None:
        return (False, "Molecule does not contain a fused steroid core (expected a subset of 4 fused carbocyclic rings: three six-membered and one five-membered).")
    
    # Within the detected steroid core, search for a candidate 17beta-hydroxy group.
    # We require a carbon in the steroid core (and specifically in the five-membered ring)
    # that has defined chirality and has a single-bonded oxygen neighbor.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:
            continue
        if atom.GetIdx() not in steroid_core_atoms:
            continue
        # Check that chirality is specified.
        if atom.GetChiralTag() == Chem.rdchem.ChiralType.CHI_UNSPECIFIED:
            continue
        # Must belong to the unique five-membered ring.
        if five_ring_indices is None or atom.GetIdx() not in five_ring_indices:
            continue
        
        # Look for -OH group: find an oxygen neighbor connected by a SINGLE bond.
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 8:
                bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
                if bond and bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                    return True, "Molecule contains a fused steroid nucleus with a candidate 17beta-hydroxy group."
                    
    return False, "No candidate 17beta-hydroxy group was detected in the steroid nucleus."

# Example usage (for testing purposes):
if __name__ == '__main__':
    # Test with one example (17beta-estradiol)
    test_smiles = "[H][C@]12CC[C@]3(C)[C@@H](O)CC[C@@]3([H])[C@]1([H])CCc1cc(O)ccc21"
    result, reason = is_17beta_hydroxy_steroid(test_smiles)
    print("Result:", result)
    print("Reason:", reason)