"""
Classifies: CHEBI:35343 17beta-hydroxy steroid
"""
"""
Classifies: CHEBI:17beta-hydroxy steroid
Definition: A 17-hydroxy steroid in which the hydroxy group at position 17 has a beta-configuration.
Heuristic update:
  1. Identify the steroid nucleus by considering only carbocyclic rings (rings with all carbons)
     of size 5 or 6.
  2. Build a “fused” ring component from those rings (rings sharing at least 2 atoms).
  3. Accept a component only if it comprises exactly 4 rings with three 6-membered and one 5-membered ring.
  4. Within that fused steroid core, search for a chiral carbon that bears an –OH (via a single bond)
     and that belongs to the sole 5-membered ring.
     
This should help reduce false positives coming from other fused cyclic systems.
"""

from rdkit import Chem

def is_17beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 17beta-hydroxy steroid based on its SMILES string.
    The algorithm first tries to detect a fused steroid nucleus (the classic four-ring system
    composed of three six-membered rings and one five-membered ring, all carbocyclic) and then 
    checks if within the five-membered ring there is a chiral carbon bearing an –OH group.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a 17beta-hydroxy steroid, False otherwise.
        str: The explanation for the classification decision.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get ring information and focus on rings of size 5 or 6 that are *carbocyclic* (only carbon atoms)
    ring_info = mol.GetRingInfo().AtomRings()
    carbocyclic_rings = []
    for ring in ring_info:
        if len(ring) in (5, 6):
            # Check that every atom in the ring is carbon (atomic number 6)
            if all(mol.GetAtomWithIdx(idx).GetAtomicNum() == 6 for idx in ring):
                carbocyclic_rings.append(set(ring))
    if not carbocyclic_rings:
        return False, "No 5- or 6-membered carbocyclic rings found; unlikely to be a steroid"

    # Build a graph among rings: each ring is a node; an edge exists if rings share at least 2 atoms.
    n = len(carbocyclic_rings)
    ring_graph = {i: set() for i in range(n)}
    for i in range(n):
        for j in range(i+1, n):
            if len(carbocyclic_rings[i].intersection(carbocyclic_rings[j])) >= 2:
                ring_graph[i].add(j)
                ring_graph[j].add(i)
    
    # Find connected components of rings in the graph.
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
    
    # Look through each connected component to see if any qualifies as a steroid core:
    # classical steroid has four fused rings: three six-membered and one five-membered ring.
    steroid_core_atoms = None
    five_ring_indices = None   # will hold set(s) corresponding to the 5-membered ring in core
    for comp in components:
        comp_rings = [carbocyclic_rings[i] for i in comp]
        # Accept only if exactly four rings are present
        if len(comp_rings) != 4:
            continue
        count_5 = sum(1 for ring in comp_rings if len(ring) == 5)
        count_6 = sum(1 for ring in comp_rings if len(ring) == 6)
        if count_5 == 1 and count_6 == 3:
            # Merge all atom indices in the component to represent the fused steroid core.
            core_atoms = set().union(*comp_rings)
            steroid_core_atoms = core_atoms
            # Also record the indices of the five-membered ring (should be only one)
            for ring in comp_rings:
                if len(ring) == 5:
                    five_ring_indices = ring
                    break
            # Use the first component matching our criteria.
            break

    if steroid_core_atoms is None:
        return False, ("Molecule does not contain a fused steroid core "
                       "(expected exactly four fused carbocyclic rings: three six-membered and one five-membered)")
    
    # Within the detected steroid core, look for a candidate 17beta-hydroxy group.
    # We require a carbon in the steroid core that:
    #   - has defined chirality (i.e. not CHI_UNSPECIFIED)
    #   - is attached by a SINGLE bond to an oxygen (–OH)
    #   - lies in the five-membered ring (assigned to be the D ring, where C17 is located)
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:
            continue  # skip non-carbons
        if atom.GetIdx() not in steroid_core_atoms:
            continue  # must be in the steroid nucleus
        # Check that chirality is specified (assume if specified then configuration is defined)
        if atom.GetChiralTag() == Chem.rdchem.ChiralType.CHI_UNSPECIFIED:
            continue
        
        # Ensure the atom is part of the five-membered ring (the candidate D ring)
        if five_ring_indices is None or atom.GetIdx() not in five_ring_indices:
            continue
        
        # Look for an -OH group: a single bonded oxygen neighbor.
        has_OH = False
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 8:
                # Check that the oxygen is connected by a single bond.
                bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
                if bond is not None and bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                    has_OH = True
                    break
        if has_OH:
            return True, "Molecule contains a fused steroid nucleus with a candidate 17beta-hydroxy group."
    
    return False, "No candidate 17beta-hydroxy group was detected in the steroid nucleus."

# Example usage (for testing only):
if __name__ == '__main__':
    # You might test with one known example such as 17beta-estradiol:
    test_smiles = "[H][C@]12CC[C@]3(C)[C@@H](O)CC[C@@]3([H])[C@]1([H])CCc1cc(O)ccc21"
    result, reason = is_17beta_hydroxy_steroid(test_smiles)
    print("Result:", result)
    print("Reason:", reason)