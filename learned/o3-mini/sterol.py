"""
Classifies: CHEBI:15889 sterol
"""
#!/usr/bin/env python
"""
Classifies: Sterol
Definition: Any 3-hydroxy steroid whose skeleton is closely related to cholestan-3-ol.
Improved strategy:
  1. Identify fused ring systems via RDKit ring info.
  2. Within each fused system, look for a candidate nucleus defined as a set of 4 rings 
     whose sizes (sorted) equal [5,6,6,6] and whose union contains ~16-19 carbon atoms.
  3. Verify that at least one of these nucleus carbons bears an external hydroxyl group.
If found, returns True and explanation; otherwise returns False and a reason.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
import itertools

def is_sterol(smiles: str):
    """
    Determines if a molecule is a sterol (a 3-hydroxy steroid) based on its SMILES string.
    
    The algorithm:
       1. Parses the SMILES and adds explicit hydrogens.
       2. Obtains the ring information from RDKit.
       3. Groups ring atoms that are fused (connected by ring bonds).
       4. In each fused component, iterates over all combinations of 4 rings.
          A valid steroid nucleus candidate must have rings whose sizes (sorted) equal [5,6,6,6] 
          and its union contains 16-19 carbon atoms.
       5. Verifies that at least one carbon in the candidate nucleus has an external hydroxyl 
          (an oxygen with at least one hydrogen attached that is not part of the nucleus).
    
    Args:
        smiles (str): The SMILES string for the molecule.
    
    Returns:
        bool: True if the molecule is classified as a sterol, False otherwise.
        str: Explanation for the classification decision.
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Adding explicit hydrogens helps in recognizing hydroxyl groups with hydrogens.
    mol = Chem.AddHs(mol)

    ring_info = mol.GetRingInfo()
    all_rings = ring_info.AtomRings()
    if not all_rings:
        return False, "No ring system present; cannot be a sterol"

    # Get set of all atoms that are parts of rings.
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
    
    # Find fused (connected) ring components.
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
    
    # For each fused ring component, try to find a steroid nucleus candidate.
    candidate_found = False
    candidate_reason = ""
    for comp in fused_components:
        # Get rings completely contained in the fused component.
        comp_rings = [r for r in all_rings if set(r).issubset(comp)]
        if len(comp_rings) < 4:
            continue  # need at least 4 rings for steroid nucleus
        
        # Try each combination of 4 rings.
        for four_rings in itertools.combinations(comp_rings, 4):
            sizes = sorted([len(r) for r in four_rings])
            # In a typical steroid nucleus, the four rings are one five-membered ring and three six-membered rings.
            if sizes != [5, 6, 6, 6]:
                continue
            
            # The union of these rings is our candidate nucleus.
            nucleus = set().union(*four_rings)
            # Count carbon atoms in this candidate nucleus.
            num_carbons = sum(1 for idx in nucleus if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
            # In a classical steroid nucleus (cholestane core), expect ~17 carbons with a small tolerance.
            if not (16 <= num_carbons <= 19):
                continue
            
            # Check for an external hydroxyl (oxygen with at least one hydrogen) attached to a nucleus carbon.
            hydroxyl_found = False
            for atom_idx in nucleus:
                atom = mol.GetAtomWithIdx(atom_idx)
                if atom.GetAtomicNum() != 6:
                    continue  # consider only carbons in the nucleus
                for nbr in atom.GetNeighbors():
                    # We require the neighbor to be oxygen and not in the nucleus.
                    if nbr.GetAtomicNum() == 8 and nbr.GetIdx() not in nucleus:
                        # Check if the oxygen is part of an –OH group (has at least one hydrogen neighbor)
                        if any(n.GetAtomicNum() == 1 for n in nbr.GetNeighbors()):
                            hydroxyl_found = True
                            break
                if hydroxyl_found:
                    break
            
            if hydroxyl_found:
                candidate_found = True
                candidate_reason = ("Found fused steroid nucleus (4 rings with sizes [5,6,6,6] and "
                                    f"approximately {num_carbons} carbons) with an externally attached hydroxyl")
                break  # nucleus candidate found in this component
            else:
                candidate_reason = ("Fused ring system candidate found with four rings ([5,6,6,6] and ~"
                                    f"{num_carbons} carbons) but missing an attached hydroxyl group")
        if candidate_found:
            break

    if candidate_found:
        return True, candidate_reason
    else:
        # If we ran through all fused components without success, provide the most-relevant reason.
        if candidate_reason:
            return False, candidate_reason
        else:
            return False, "No steroid nucleus identified with the required ring pattern"

# Example usage
if __name__ == "__main__":
    # Test with one known sterol: (25R)-3beta,4beta-dihydroxycholest-5-en-26-oic acid
    test_smiles = "C[C@H](CCC[C@@H](C)C(O)=O)[C@H]1CC[C@H]2[C@@H]3CC=C4[C@@H](O)[C@@H](O)CC[C@]4(C)[C@H]3CC[C@]12C"
    res, reason = is_sterol(test_smiles)
    print("Is sterol:", res)
    print("Reason:", reason)