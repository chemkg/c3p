"""
Classifies: CHEBI:36835 3alpha-hydroxy steroid
"""
"""
Classifies: 3alpha-hydroxy steroid
Definition: A 3-hydroxy steroid in which the 3-hydroxy substituent is in the alpha-position.
This improved version first identifies a fused ring system that is (a) large (at least 3 rings),
(b) contains at least one 5-membered ring and mostly 5- or 6-membered rings, and (c) is almost
exclusively composed of carbon centers – as expected for a cyclopentanoperhydrophenanthrene steroid nucleus.
Then it checks that at least one alpha-OH (indicated by the [C@@H](O) SMARTS) is attached to one
of those fused ring atoms.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_3alpha_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 3alpha-hydroxy steroid based on its SMILES string.
    Improved heuristics:
      1. Identify a fused ring system from ring information. The steroid nucleus is known to be a
         single fused set of rings (typically four rings) that is largely all-carbon. Here we extract
         the largest connected set (by ring fusion, where two rings are deemed fused if they share at least 2 atoms).
      2. Require that the fused component contains at least three rings and at least one five-membered ring.
      3. Check that most atoms in the fused ring system are carbons (since steroid backbones are carboskeletal).
      4. Look for at least one [C@@H](O) substructure (the alpha hydroxyl proxy) whose carbon atom is part
         of the fused ring nucleus.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule is classified as a 3alpha-hydroxy steroid, else False.
        str: Reason for the classification decision.
    """
    # Parse the SMILES and assign stereochemistry
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    
    # Get ring information from the molecule.
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    if not atom_rings:
        return False, "No rings found, not a steroid"
    
    # Build a graph of rings: each ring is a node; add an edge between rings sharing >= 2 atoms.
    n_rings = len(atom_rings)
    # Create an adjacency list for ring indices
    adj = {i: set() for i in range(n_rings)}
    for i in range(n_rings):
        for j in range(i+1, n_rings):
            if len(set(atom_rings[i]).intersection(atom_rings[j])) >= 2:
                adj[i].add(j)
                adj[j].add(i)
    
    # Find connected components (fused ring sets) using DFS.
    seen = set()
    fused_components = []
    for i in range(n_rings):
        if i in seen:
            continue
        stack = [i]
        component = set()
        while stack:
            node = stack.pop()
            if node in component:
                continue
            component.add(node)
            for neighbor in adj[node]:
                if neighbor not in component:
                    stack.append(neighbor)
        seen.update(component)
        fused_components.append(component)
    
    # Choose the largest fused set (by number of rings)
    largest_component = max(fused_components, key=lambda comp: len(comp))
    # Build a set of all atom indices in this fused ring system
    fused_atoms = set()
    for idx in largest_component:
        fused_atoms.update(atom_rings[idx])
    
    # Check that the fused system has at least 3 rings.
    if len(largest_component) < 3:
        return False, "Fused ring system too small to be a steroid nucleus"
    
    # For steroid nucleus, we expect mostly rings that are 5- or 6-membered.
    num_5membered = sum(1 for i in largest_component if len(atom_rings[i]) == 5)
    num_6membered = sum(1 for i in largest_component if len(atom_rings[i]) == 6)
    if num_5membered < 1:
        return False, "Fused ring system missing a five-membered ring typical of a steroid nucleus"
    # Optionally: require that all rings in the component have 5 or 6 atoms
    for i in largest_component:
        if len(atom_rings[i]) not in (5, 6):
            return False, "Fused ring system contains rings not typical for a steroid nucleus (only 5- or 6-membered allowed)"
    
    # Check that the fused nucleus is mostly carbon (at least 90% of atoms should be carbon).
    nucleus_atom_count = len(fused_atoms)
    carbon_count = 0
    for atom_idx in fused_atoms:
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom.GetAtomicNum() == 6:
            carbon_count += 1
    if nucleus_atom_count == 0 or (carbon_count / nucleus_atom_count) < 0.9:
        return False, "Fused ring nucleus is not predominantly carbon as expected for a steroid"
    
    # Now, look for an alpha hydroxyl group.
    # We use the SMARTS "[C@@H](O)" as a proxy for an alpha-OH.
    alpha_oh_smarts = "[C@@H](O)"
    alpha_oh_pattern = Chem.MolFromSmarts(alpha_oh_smarts)
    if alpha_oh_pattern is None:
        return False, "Invalid SMARTS for alpha hydroxyl group"
    
    matches = mol.GetSubstructMatches(alpha_oh_pattern)
    if not matches:
        return False, "No [C@@H](O) substructure found (alpha hydroxyl group not detected)"
    
    # Check that at least one [C@@H](O) match occurs on an atom that is inside the fused nucleus.
    alpha_on_nucleus = False
    for match in matches:
        # match[0] is the carbon bearing the OH.
        if match[0] in fused_atoms:
            alpha_on_nucleus = True
            break
    if not alpha_on_nucleus:
        return False, "Alpha hydroxyl group does not reside on the fused steroid nucleus"
    
    return True, "Molecule has a fused (5,6)-ring nucleus with predominantly carbon atoms and a 3α-hydroxy substituent"

# Testing block:
if __name__ == "__main__":
    # Example: 5b-Cholestane-3a,7a,12a,23-tetrol (true positive)
    test_smiles = "C1C[C@H](C[C@]2(C[C@H]([C@@]3([C@](C[C@@H]([C@]4([C@]3(CC[C@@]4([C@@](CC(CC(C)C)O)(C)[H])[H])[H])C)O)([C@@]12C)[H])[H])O)[H])O"
    result, reason = is_3alpha_hydroxy_steroid(test_smiles)
    print("Result:", result)
    print("Reason:", reason)