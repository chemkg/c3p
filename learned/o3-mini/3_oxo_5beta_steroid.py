"""
Classifies: CHEBI:1624 3-oxo-5beta-steroid
"""
"""
Classifies: 3-oxo-5beta-steroid, defined as 'Any 3-oxo steroid that has beta- configuration at position 5.' 
Heuristic criteria revised:
  1. Identify all ring systems from the input molecule.
  2. Build a graph of fused rings (rings sharing at least 2 atoms) and search for a connected component
     that has exactly 4 rings and whose ring sizes are exactly three six-membered rings and one five-membered ring.
     (This represents the canonical steroid nucleus, i.e. cyclopentanoperhydrophenanthrene.)
  3. Check that at least one of the carbons in that candidate nucleus bears a double-bonded oxygen (a ketone) 
     – which should be a ring-bound carbonyl (the "3-oxo" feature).
  4. As a proxy for the beta configuration at carbon-5, require that the SMILES string includes at least one "@@" annotation.
     
Note: Since steroid numbering is nontrivial from a SMILES string, these conditions are all heuristic.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.rdchem import BondType

def is_3_oxo_5beta_steroid(smiles: str):
    """
    Determines if a molecule appears to be a 3-oxo-5beta-steroid based on its SMILES string.
    
    The heuristic:
      (a) Parse the molecule and extract all ring systems.
      (b) Build a connectivity graph among rings that share at least 2 atoms. Then search for a connected
          set of rings that has exactly 4 fused rings, consisting of 3 six-membered rings and 1 five-membered ring.
      (c) Check that within the candidate steroid nucleus there is at least one ring-bound ketone (C=O).
      (d) Verify that the original SMILES string contains an explicit chiral tag ("@@") that is assumed to
          indicate the beta configuration at position 5.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule fits the heuristic criteria for a 3-oxo-5beta-steroid, False otherwise.
        str: Explanation for the decision.
    """
    # Parse the SMILES string into a molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get ring information: each ring is returned as a tuple of atom indices.
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    if not rings:
        return False, "No rings found in the molecule"
    
    # Build a connectivity graph of rings.
    n_rings = len(rings)
    ring_nodes = list(range(n_rings))
    ring_sets = [set(r) for r in rings]  # convert tuples to sets for intersection
    connectivity = {i: set() for i in ring_nodes}
    for i in ring_nodes:
        for j in ring_nodes:
            if i < j and len(ring_sets[i].intersection(ring_sets[j])) >= 2:
                connectivity[i].add(j)
                connectivity[j].add(i)
    
    # Find connected components in the ring graph.
    visited = set()
    components = []
    for node in ring_nodes:
        if node not in visited:
            stack = [node]
            comp = set()
            while stack:
                current = stack.pop()
                if current not in visited:
                    visited.add(current)
                    comp.add(current)
                    stack.extend(list(connectivity[current] - visited))
            components.append(comp)
    
    # Look for a connected component that represents the steroid nucleus:
    # It must have exactly 4 rings, with exactly 3 six-membered rings and 1 five-membered ring.
    candidate_component = None
    for comp in components:
        if len(comp) == 4:
            six_count = 0
            five_count = 0
            for idx in comp:
                ring_size = len(rings[idx])
                if ring_size == 6:
                    six_count += 1
                elif ring_size == 5:
                    five_count += 1
            if six_count == 3 and five_count == 1:
                candidate_component = comp
                break
    if candidate_component is None:
        return False, ("Fused ring system does not match typical steroid nucleus. "
                       "Expected exactly 4 fused rings (3 six-membered and 1 five-membered) but none found.")
    
    # Create the set of atoms that are part of the candidate steroid nucleus.
    candidate_atoms = set()
    for idx in candidate_component:
        candidate_atoms.update(rings[idx])
    
    # Search for a ketone group (C=O) where the carbon is in the candidate nucleus.
    ketone_found = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6 and atom.GetIdx() in candidate_atoms:
            # Check if this carbon is double-bonded to an oxygen.
            for bond in atom.GetBonds():
                # Check the bond type using RDKit's BondType enum.
                if bond.GetBondType() == BondType.DOUBLE:
                    nbr = bond.GetOtherAtom(atom)
                    if nbr.GetAtomicNum() == 8:
                        ketone_found = True
                        break
            if ketone_found:
                break
    if not ketone_found:
        return False, "No ring-bound ketone group found within the candidate steroid nucleus (missing 3-oxo feature)"
    
    # Check for chiral annotation indicating beta configuration.
    # We require that the original SMILES contains the '@@' symbol.
    if "@@" not in smiles:
        return False, "No '@@' chiral annotation found in SMILES (missing beta configuration indicator)"
    
    return True, "Molecule appears to be a 3-oxo steroid with beta configuration based on heuristic criteria"

# Example usage: (uncomment the following lines to test)
if __name__ == "__main__":
    # A test SMILES for 17beta-hydroxy-5beta-estran-3-one
    test_smiles = "C[C@]12CC[C@H]3[C@@H](CC[C@@H]4CC(=O)CC[C@H]34)[C@@H]1CC[C@@H]2O"
    result, reason = is_3_oxo_5beta_steroid(test_smiles)
    print(result, reason)