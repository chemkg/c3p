"""
Classifies: CHEBI:35343 17beta-hydroxy steroid
"""
#!/usr/bin/env python
"""
Classifies: 17beta-hydroxy steroid
Definition: A 17-hydroxy steroid in which the hydroxy group at position 17 has a beta-configuration.
Heuristic:
  1. Identify the fused steroid nucleus by extracting rings of size 5 and 6 and looking for a connected cluster
     having at least one 5-membered and three 6-membered rings.
  2. Within that fused core, search for a chiral carbon that bears an –OH and that belongs to a 5-membered ring.
     (We assume that if a stereocenter is defined on such a candidate it is interpreted as beta.)
"""

from rdkit import Chem

def is_17beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 17beta-hydroxy steroid based on its SMILES string.
    The algorithm checks that the molecule contains a fused steroid nucleus (at least one 5-membered and three 6-membered rings
    that are interconnected) and that within that nucleus there is a chiral carbon in a five-membered ring bearing an –OH group.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule is a 17beta-hydroxy steroid, False otherwise.
        str: The explanation for the classification decision.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Obtain ring information (list of tuples of atom indices)
    rings = mol.GetRingInfo().AtomRings()
    
    # Focus only on rings of size 5 or 6 (typical for steroid nuclei)
    steroid_rings = [set(ring) for ring in rings if len(ring) in (5, 6)]
    if not steroid_rings:
        return False, "No 5- or 6-membered rings found; unlikely to be a steroid"
    
    # Build a graph of rings: nodes = each ring; an edge between two rings exists if they share at least 2 atoms.
    n = len(steroid_rings)
    ring_graph = {i: set() for i in range(n)}
    for i in range(n):
        for j in range(i+1, n):
            if len(steroid_rings[i].intersection(steroid_rings[j])) >= 2:
                ring_graph[i].add(j)
                ring_graph[j].add(i)
    
    # Find connected components of rings from the graph.
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
    
    # Identify a component that approximates the fused steroid core:
    # We want at least one 5-membered ring and at least three 6-membered rings in the connected set.
    steroid_component_atoms = None
    for comp in components:
        comp_rings = [steroid_rings[i] for i in comp]
        count_5 = sum(1 for ring in comp_rings if len(ring) == 5)
        count_6 = sum(1 for ring in comp_rings if len(ring) == 6)
        if count_5 >= 1 and count_6 >= 3:
            # Merge all atom indices in this component to represent the steroid core.
            steroid_component_atoms = set().union(*comp_rings)
            break

    if steroid_component_atoms is None:
        return False, ("Molecule does not contain a fused steroid core "
                       "(expected at least one 5-membered and three 6-membered fused rings)")
    
    # Now search within the steroid core for a candidate 17beta-hydroxy group:
    # Look for a chiral carbon (with explicit stereochemistry) that is in the steroid core,
    # has an -OH (an oxygen via single bond) and belongs to at least one 5-membered ring.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:
            continue
        # Require that the atom has defined chirality
        if atom.GetChiralTag() == Chem.rdchem.ChiralType.CHI_UNSPECIFIED:
            continue
        # Ensure the atom is part of the steroid core
        if atom.GetIdx() not in steroid_component_atoms:
            continue
        # Check for an attached hydroxyl group (oxygen connected by a SINGLE bond)
        has_OH = False
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 8:
                bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
                if bond is not None and bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                    has_OH = True
                    break
        if not has_OH:
            continue
        # Check that this candidate atom is in at least one 5-membered ring of the molecule
        in_five_ring = False
        for ring in rings:
            if len(ring) == 5 and atom.GetIdx() in ring:
                # Optionally, make sure the entire ring is part of the steroid core
                if steroid_component_atoms.issuperset(ring):
                    in_five_ring = True
                    break
        if in_five_ring:
            return True, "Molecule contains a fused steroid nucleus with a candidate 17beta-hydroxy group."
    
    return False, "No candidate 17beta-hydroxy group was detected in the steroid nucleus."

# Example usage:
if __name__ == '__main__':
    # Test with one example: 17beta-estradiol 
    test_smiles = "[H][C@]12CC[C@]3(C)[C@@H](O)CC[C@@]3([H])[C@]1([H])CCc1cc(O)ccc21"
    result, reason = is_17beta_hydroxy_steroid(test_smiles)
    print("Result:", result)
    print("Reason:", reason)