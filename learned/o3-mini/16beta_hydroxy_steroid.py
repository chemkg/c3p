"""
Classifies: CHEBI:17354 16beta-hydroxy steroid
"""
#!/usr/bin/env python
"""
Classifies: 16beta-hydroxy steroid
A 16-beta-hydroxy steroid is defined as a steroid (with a fused tetracyclic nucleus typical of steroids—
namely three six-membered rings and one five-membered ring, normally called the A–B–C–D rings)
in which an –OH group is attached at the position (presumably C16 on the D ring) and drawn with explicit stereochemistry (assumed beta).
This heuristic implementation:
  1. Parses the molecule and adds explicit hydrogens.
  2. Retrieves all rings and builds a graph of rings that are "fused" (i.e. share ≥2 atoms).
  3. Searches for a connected component (fused ring system) that contains at least 1 five-membered ring and 3 six-membered rings.
  4. Within such a component, looks for a carbon (an atom in a five-membered ring) that has an –OH neighbor (an oxygen bonded to at least one H)
     and whose chirality is explicitly defined (assumed to indicate beta configuration).
If these checks are met, the molecule is classified as a 16beta-hydroxy steroid.
Otherwise, the function returns False with a reason.
Note: This is a heuristic approximate method.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_16beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 16beta-hydroxy steroid based on its SMILES string.
    That is, it contains a fused steroid nucleus (typically 3 six-membered rings and 1 five-membered ring fused)
    and an -OH group attached to a carbon in the five-membered ring with defined stereochemistry.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule is a 16beta-hydroxy steroid, False otherwise.
        str: Reason for classification.
    """
    # Parse molecule from SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens to see –OH groups clearly and preserve stereochemistry.
    mol = Chem.AddHs(mol)
    
    # Obtain ring information.
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    if not rings:
        return False, "No rings found in molecule; not a steroid nucleus"
    
    # Build list of rings with their sizes and set versions of atom indices.
    rings_list = []
    for ring in rings:
        rings_list.append(set(ring))
    
    # Build a graph of rings: nodes indices (0,1,2,...) and add an edge if two rings share at least 2 atoms.
    n_rings = len(rings_list)
    ring_graph = {i: set() for i in range(n_rings)}
    for i in range(n_rings):
        for j in range(i+1, n_rings):
            if len(rings_list[i].intersection(rings_list[j])) >= 2:
                ring_graph[i].add(j)
                ring_graph[j].add(i)
    
    # Find connected components (fused ring systems) in the ring graph.
    visited = set()
    components = []
    for i in range(n_rings):
        if i not in visited:
            # Use DFS to collect component nodes.
            stack = [i]
            comp = set()
            while stack:
                node = stack.pop()
                if node not in comp:
                    comp.add(node)
                    for nei in ring_graph[node]:
                        if nei not in comp:
                            stack.append(nei)
            visited |= comp
            components.append(comp)
    
    # Look for a component that is steroid-like i.e. it contains at least one 5-membered ring and at least three 6-membered rings.
    steroid_component = None
    for comp in components:
        count_5 = 0
        count_6 = 0
        for idx in comp:
            size = len(rings_list[idx])
            if size == 5:
                count_5 += 1
            elif size == 6:
                count_6 += 1
        if count_5 >= 1 and count_6 >= 3:
            steroid_component = comp
            break
    
    if steroid_component is None:
        return False, "Molecule does not have a fused ring system with the typical steroid nucleus (needs at least 1 five-membered and 3 six-membered fused rings)"
    
    # Get the union of atom indices that are part of the fused steroid nucleus.
    steroid_atoms = set()
    five_membered_rings = []  # focus on candidate 5-membered rings from the steroid nucleus
    for idx in steroid_component:
        ring_atoms = rings_list[idx]
        steroid_atoms |= ring_atoms
        if len(ring_atoms) == 5:
            five_membered_rings.append(ring_atoms)
    
    # Now scan atoms in the fused steroid nucleus that belong to one of the five-membered rings.
    candidate_found = False
    candidate_has_stereo = False
    candidate_reason = ""
    for ring in five_membered_rings:
        for atom_idx in ring:
            atom = mol.GetAtomWithIdx(atom_idx)
            # We focus on carbon atoms.
            if atom.GetAtomicNum() != 6:
                continue
            # Look for neighbor oxygen that is –OH.
            for nbr in atom.GetNeighbors():
                if nbr.GetAtomicNum() == 8:
                    # Confirm that the oxygen is part of a hydroxyl (has at least one H neighbor)
                    if any(h.GetAtomicNum() == 1 for h in nbr.GetNeighbors()):
                        candidate_found = True
                        # Check chirality on the carbon. A defined chiral tag is taken as a proxy for beta configuration.
                        chiral = atom.GetChiralTag()
                        if chiral != Chem.rdchem.ChiralType.CHI_UNSPECIFIED:
                            candidate_has_stereo = True
                            candidate_reason = ("Found fused steroid nucleus (at least 3 six-membered and 1 five-membered fused rings) "
                                                "with an -OH group attached to a carbon in a five-membered ring that has defined stereochemistry "
                                                "(assumed 16beta-hydroxy).")
                            break
                        else:
                            candidate_reason = ("Found fused steroid nucleus with an -OH group on a carbon in a five-membered ring, "
                                                "but the carbon lacks explicit stereochemistry; cannot confirm beta-configuration.")
            if candidate_found and candidate_has_stereo:
                break
        if candidate_found and candidate_has_stereo:
            break

    if candidate_found:
        if candidate_has_stereo:
            return True, candidate_reason
        else:
            return False, candidate_reason
    else:
        return False, "Could not find an -OH group on a carbon (with defined stereochemistry) in any five-membered ring of the fused steroid nucleus (position 16 candidate not found)"

# For testing purposes:
if __name__ == "__main__":
    # Example SMILES strings
    test_smiles = [
        "C[C@]12CC[C@H]3[C@@H](CCC4=CC(=O)CC[C@]34C)[C@@H]1C[C@H](O)[C@@H]2O",  # 16beta-hydroxytestosterone
        "C1=C2C(CC[C@]3([C@@]4(C[C@@H]([C@@H]([C@]4(CC[C@@]32[H])C)O)O)[H])[H])=CC(=C1)O",  # 16beta-hydroxyestradiol (example)
    ]
    for s in test_smiles:
        result, reason = is_16beta_hydroxy_steroid(s)
        print("SMILES:", s)
        print("Classification:", result)
        print("Reason:", reason)
        print("-----")