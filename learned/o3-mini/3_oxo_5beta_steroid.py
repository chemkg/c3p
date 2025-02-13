"""
Classifies: CHEBI:1624 3-oxo-5beta-steroid
"""
"""
Classifies: 3-oxo-5beta-steroid, defined as 'Any 3-oxo steroid that has beta- configuration at position 5.'
Heuristic criteria:
  1. The molecule must have a steroid-like fused ring system. We identify rings via RDKit’s ring info,
     then build a graph of rings that share at least 2 atoms. The largest connected component is assumed to be the steroid nucleus.
     A typical steroid core is composed of 4 fused rings (three six-membered rings and one five-membered ring).
  2. The molecule should contain a ketone group (C=O) in which the carbon is inside a ring (to provide the “3-oxo” feature).
  3. The input SMILES should include the “@@” stereochemical marker as a proxy for the beta configuration at position 5.
Note: These checks are heuristic and might produce false positives or negatives.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_3_oxo_5beta_steroid(smiles: str):
    """
    Determines if a molecule belongs to the 3-oxo-5beta-steroid class based on its SMILES string.
    
    The function performs a series of heuristic checks:
      (a) Parses the molecule and gets its ring information.
      (b) Constructs a graph of rings that are fused (sharing at least 2 atoms)
          and verifies that the largest connected component has at least four rings,
          including at least three six-membered rings and one five-membered ring.
      (c) Searches for a ketone (C=O) where the carbon is part of a ring.
      (d) Looks in the original SMILES string for an explicit "@@" chiral annotation.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule passes the heuristic checks, False otherwise.
        str: Explanation for the decision.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get ring information from the molecule
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()  # Each ring is a tuple of atom indices.
    if not rings:
        return False, "No rings found in the molecule"
    
    # Build a "ring connectivity" graph: each ring is a node;
    # connect two rings if they share at least two atoms.
    ring_nodes = list(range(len(rings)))
    # Represent rings as sets for easier intersection
    ring_sets = [set(r) for r in rings]
    connectivity = {i: set() for i in ring_nodes}
    for i in ring_nodes:
        for j in ring_nodes:
            if i < j:
                if len(ring_sets[i].intersection(ring_sets[j])) >= 2:
                    connectivity[i].add(j)
                    connectivity[j].add(i)
    
    # Find connected components in the ring graph
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
                    stack.extend(connectivity[current] - visited)
            components.append(comp)
    
    # Identify the largest connected component (candidate steroid nucleus)
    largest_comp = max(components, key=lambda comp: len(comp))
    # Count the rings by size in that component.
    six_count = 0
    five_count = 0
    for idx in largest_comp:
        ring_size = len(rings[idx])
        if ring_size == 6:
            six_count += 1
        elif ring_size == 5:
            five_count += 1
    # Typical steroid nucleus (cyclopentanoperhydrophenanthrene) has 3 six-membered and 1 five-membered ring.
    if len(largest_comp) < 4:
        return False, ("Fused ring system too small. "
                       "Expected at least 4 contiguous rings as a steroid core, found {}.".format(len(largest_comp)))
    if six_count < 3 or five_count < 1:
        return False, ("Ring sizes in fused system do not match steroid pattern "
                       "(requires at least three 6-membered and one 5-membered rings; found {} six-membered and {} five-membered).".format(six_count, five_count))
    
    # Check for a ketone group (C=O) where the carbon is in a ring:
    # Loop over all atoms to find a carbon with a double-bonded oxygen and verify the carbon is in a ring.
    ketone_found = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:  # carbon
            # Check if atom is in a ring
            if mol.GetRingInfo().NumAtomRings(atom.GetIdx()) > 0:
                for bond in atom.GetBonds():
                    # Look for a double bond to oxygen
                    nbr = bond.GetOtherAtom(atom)
                    if nbr.GetAtomicNum() == 8 and bond.GetBondTypeAsDouble() == 2:
                        ketone_found = True
                        break
            if ketone_found:
                break
    if not ketone_found:
        return False, "No ring-bound ketone group found (missing 3-oxo feature)"
    
    # Check for chiral annotation for beta configuration.
    # We use the presence of "@@" in the SMILES string as a proxy.
    if "@@" not in smiles:
        return False, "No '@@' stereochemical annotation found (missing beta configuration indicator)"
    
    return True, "Molecule appears to be a 3-oxo steroid with beta configuration based on heuristic criteria"

# Example usage: (uncomment the following lines to test)
if __name__ == "__main__":
    # A test SMILES for 17beta-hydroxy-5beta-estran-3-one
    test_smiles = "C[C@]12CC[C@H]3[C@@H](CC[C@@H]4CC(=O)CC[C@H]34)[C@@H]1CC[C@@H]2O"
    result, reason = is_3_oxo_5beta_steroid(test_smiles)
    print(result, reason)