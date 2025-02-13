"""
Classifies: CHEBI:15889 sterol
"""
#!/usr/bin/env python
"""
Classifies: Sterol
Definition: Any 3-hydroxy steroid whose skeleton is closely related to cholestan-3-ol.
This implementation:
  1. Parses the SMILES and adds explicit hydrogens.
  2. Retrieves the ring information from RDKit.
  3. Groups ring atoms into fused components.
  4. For each fused component, iterates over all combinations of 4 rings.
     For a candidate steroid nucleus the four rings should have sizes [5,6,6,6] when sorted,
     the union (candidate nucleus) should have between 16 and 19 carbon atoms, and—importantly—
     the four rings must be fused in a steroid‐like pattern. This connectivity is verified by
     requiring that exactly three pairs of rings (out of six possible pairings) share at least
     2 atoms and that the resulting graph (with rings as nodes and fused pairs as edges) is connected.
  5. Finally, at least one carbon in the candidate steroid nucleus must bear an external –OH group.
If all these conditions are met, returns True with explanation; otherwise False.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
import itertools

def is_sterol(smiles: str):
    """
    Determines if a molecule is a sterol (a 3-hydroxy steroid) based on its SMILES string.
    
    The algorithm:
       1. Parse the SMILES, and add explicit hydrogens.
       2. Obtain ring information and group rings into fused (connected) components using ring bonds.
       3. For each fused component, consider all combinations of four rings. For each combination:
             a. Check that the sorted ring sizes equal [5,6,6,6].
             b. Compute the union of atoms in these rings and count carbon atoms in it.
                The classical steroid nucleus is expected to have between 16 and 19 carbons.
             c. Verify that the four rings are fused in the proper steroid‐like connectivity.
                We require that exactly three pairs of rings (out of the six possible pairs) share
                at least 2 atoms and that the connectivity graph is connected (i.e. it forms a tree).
             d. Check that at least one carbon of the nucleus is attached to a hydroxyl group –
                that is, an oxygen (with at least one hydrogen) that is not part of the nucleus.
       4. If such a valid candidate is found, return True with an explanation;
          otherwise return False and the most relevant reason.
          
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if classified as a sterol, False otherwise.
        str: Explanation for the decision.
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Adding explicit hydrogens will help in detecting hydroxyl groups.
    mol = Chem.AddHs(mol)
    
    ring_info = mol.GetRingInfo()
    all_rings = ring_info.AtomRings()
    if not all_rings:
        return False, "No ring system present; cannot be a sterol"
    
    # Get a set of all atom indices that are in some ring.
    ring_atom_idxs = {atom.GetIdx() for atom in mol.GetAtoms() if atom.IsInRing()}
    
    # Build a graph among ring atoms using bonds that are in rings.
    ring_graph = {idx: set() for idx in ring_atom_idxs}
    for bond in mol.GetBonds():
        if bond.IsInRing():
            i1 = bond.GetBeginAtomIdx()
            i2 = bond.GetEndAtomIdx()
            if i1 in ring_graph and i2 in ring_graph:
                ring_graph[i1].add(i2)
                ring_graph[i2].add(i1)
                
    # Find fused (i.e. connected) groups of ring atoms.
    visited = set()
    fused_components = []
    for idx in ring_atom_idxs:
        if idx not in visited:
            comp = set()
            stack = [idx]
            while stack:
                current = stack.pop()
                if current not in visited:
                    visited.add(current)
                    comp.add(current)
                    stack.extend(ring_graph[current] - visited)
            fused_components.append(comp)
    
    candidate_found = False
    candidate_reason = ""
    
    # Go through each fused component and try to find a valid steroid nucleus candidate.
    for comp in fused_components:
        # Select rings completely contained in this fused component.
        comp_rings = [set(r) for r in all_rings if set(r).issubset(comp)]
        if len(comp_rings) < 4:
            continue  # Need at least four rings.
        
        # Iterate over all combinations of four rings.
        for four_rings in itertools.combinations(comp_rings, 4):
            sizes = sorted([len(r) for r in four_rings])
            # Check for classical steroid nucleus: one five-membered and three six-membered rings.
            if sizes != [5, 6, 6, 6]:
                continue
            
            # Connectivity check: there are 4 rings; expect that exactly 3 pairs share at least 2 atoms,
            # and the connectivity graph (rings as nodes, edge if intersection >=2) should be connected.
            edges = []
            for i, j in itertools.combinations(range(4), 2):
                inter = four_rings[i] & four_rings[j]
                if len(inter) >= 2:
                    edges.append((i, j))
            if len(edges) != 3:
                continue  # Not the steroid-like fusion pattern.
            # Build a graph from these four rings and check connectivity.
            conn_graph = {i: set() for i in range(4)}
            for i, j in edges:
                conn_graph[i].add(j)
                conn_graph[j].add(i)
            # Check connectivity by simple DFS.
            stack = [0]
            seen = set()
            while stack:
                node = stack.pop()
                if node not in seen:
                    seen.add(node)
                    stack.extend(conn_graph[node] - seen)
            if len(seen) != 4:
                continue
            
            # Now, candidate nucleus is the union of the four rings.
            nucleus = set().union(*four_rings)
            num_carbons = sum(1 for idx in nucleus if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
            if not (16 <= num_carbons <= 19):
                continue
            
            # Look for an external hydroxyl group.
            hydroxyl_found = False
            for atom_idx in nucleus:
                atom = mol.GetAtomWithIdx(atom_idx)
                if atom.GetAtomicNum() != 6:  # Only consider nucleus carbons.
                    continue
                for nbr in atom.GetNeighbors():
                    # External oxygen (not in the nucleus) with at least one hydrogen.
                    if nbr.GetAtomicNum() == 8 and nbr.GetIdx() not in nucleus:
                        if any(n.GetAtomicNum() == 1 for n in nbr.GetNeighbors()):
                            hydroxyl_found = True
                            break
                if hydroxyl_found:
                    break
            
            if hydroxyl_found:
                candidate_found = True
                candidate_reason = ("Found fused steroid nucleus (4 rings with sizes [5,6,6,6] and "
                                    f"approximately {num_carbons} carbons) with an externally attached hydroxyl")
                break  # Found a valid candidate.
            else:
                candidate_reason = ("Fused ring system candidate found with four rings ([5,6,6,6] and ~"
                                    f"{num_carbons} carbons) but missing an attached hydroxyl group")
        if candidate_found:
            break

    if candidate_found:
        return True, candidate_reason
    else:
        if candidate_reason:
            return False, candidate_reason
        else:
            return False, "No steroid nucleus identified with the required ring pattern"

# Example usage (for testing):
if __name__ == "__main__":
    # Test with one known sterol: (25R)-3beta,4beta-dihydroxycholest-5-en-26-oic acid
    test_smiles = "C[C@H](CCC[C@@H](C)C(O)=O)[C@H]1CC[C@H]2[C@@H]3CC=C4[C@@H](O)[C@@H](O)CC[C@]4(C)[C@H]3CC[C@]12C"
    result, explanation = is_sterol(test_smiles)
    print("Is sterol:", result)
    print("Reason:", explanation)