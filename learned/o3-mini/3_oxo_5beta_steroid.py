"""
Classifies: CHEBI:1624 3-oxo-5beta-steroid
"""
"""
Classifies: 3-oxo-5beta-steroid, defined as 'Any 3-oxo steroid that has beta- configuration at position 5.'

Heuristic criteria:
  1. Extract all ring systems from the input molecule.
  2. For every combination of 4 rings, check if:
       - They are fused (using edge connectivity: two rings are fused if they share at least 2 atoms)
       - They match the steroid nucleus: exactly 3 six-membered rings and 1 five-membered ring.
  3. Ensure that within the candidate nucleus, at least one carbon atom has a double bond to an oxygen (the 3-oxo feature).
  4. As a proxy for the beta configuration at position 5, require that the SMILES contains an "@@" chiral tag.
  
Since steroid numbering is nontrivial from a SMILES string, these conditions are heuristic.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.rdchem import BondType
import itertools

def is_3_oxo_5beta_steroid(smiles: str):
    """
    Determines if a molecule appears to be a 3-oxo-5beta-steroid based on its SMILES string.
    
    The heuristic:
      (a) Parse the molecule and extract its rings.
      (b) From all rings, try every combination of four rings that are mutually fused 
          (adjacent if they share at least 2 atoms) and that contain exactly 3 six-membered rings and 1 five-membered ring.
      (c) Check that within the union of these ring atoms there is at least one carbon with a double-bonded oxygen.
      (d) Verify that the original SMILES string contains an explicit '@@' annotation indicating beta configuration.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule fits the heuristic for a 3-oxo-5beta-steroid; False otherwise.
        str: Explanation for the decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Obtain ring information (each ring is a tuple of atom indices).
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    if not rings:
        return False, "No rings found in the molecule"
    
    # Helper function: Given a set of rings (list of tuples), check if the induced subgraph (based on sharing >=2 atoms)
    # is connected.
    def is_connected(ring_indices):
        # Create a connectivity dictionary for rings in the combination.
        n = len(ring_indices)
        # Map local index to the ring's global index.
        local2global = list(ring_indices)
        # Build connectivity: two rings are connected if they share at least 2 atoms.
        connectivity = {i: set() for i in range(n)}
        for i in range(n):
            for j in range(i+1, n):
                # Use set intersection on atom indices.
                if len(set(rings[local2global[i]]).intersection(rings[local2global[j]])) >= 2:
                    connectivity[i].add(j)
                    connectivity[j].add(i)
        # Now do a simple connectivity (BFS) over local indices.
        visited = set()
        stack = [0]
        while stack:
            current = stack.pop()
            if current not in visited:
                visited.add(current)
                stack.extend(connectivity[current] - visited)
        return len(visited) == n

    # Now try every combination of 4 rings.
    candidate_atoms = None
    found_candidate = False
    for combo in itertools.combinations(range(len(rings)), 4):
        # Check ring sizes: count six-membered and five-membered rings.
        six_count = sum(1 for idx in combo if len(rings[idx]) == 6)
        five_count = sum(1 for idx in combo if len(rings[idx]) == 5)
        if six_count != 3 or five_count != 1:
            continue  # does not match steroid nucleus
        # Check that these rings form a connected set based on our fusion criteria.
        if not is_connected(combo):
            continue
        # If we get here, we have a candidate steroid nucleus.
        candidate_atoms = set()
        for idx in combo:
            candidate_atoms.update(rings[idx])
        found_candidate = True
        break  # we found a valid candidate nucleus; no need to search further
    if not found_candidate:
        return False, ("Fused ring system does not match typical steroid nucleus. "
                       "Expected a set of 4 fused rings (3 six-membered and 1 five-membered) within the overall ring list.")
    
    # Check for the 3-oxo feature: a ring-bound carbonyl (C=O) where the carbon lies in the candidate nucleus.
    ketone_found = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6 and atom.GetIdx() in candidate_atoms:
            # Check bonds for a double bond to oxygen.
            for bond in atom.GetBonds():
                if bond.GetBondType() == BondType.DOUBLE:
                    nbr = bond.GetOtherAtom(atom)
                    if nbr.GetAtomicNum() == 8:
                        ketone_found = True
                        break
            if ketone_found:
                break
    if not ketone_found:
        return False, "No ring-bound ketone group (3-oxo) found within the candidate steroid nucleus"
    
    # Check for beta configuration proxy: look for the '@@' annotation in the original SMILES.
    if "@@" not in smiles:
        return False, "No '@@' chiral annotation found in SMILES (missing beta configuration indicator)"
    
    return True, "Molecule appears to be a 3-oxo steroid with beta configuration based on heuristic criteria"

# Example usage: (uncomment the following code to run test examples)
if __name__ == "__main__":
    test_smiles = "C[C@]12CC[C@H]3[C@@H](CC[C@@H]4CC(=O)CC[C@H]34)[C@@H]1CC[C@@H]2O"  # 17beta-hydroxy-5beta-estran-3-one
    result, reason = is_3_oxo_5beta_steroid(test_smiles)
    print(result, reason)