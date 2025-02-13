"""
Classifies: CHEBI:36835 3alpha-hydroxy steroid
"""
"""
Classifies: 3α-hydroxy steroid
Definition: A 3–hydroxy steroid in which the 3–hydroxy substituent is in the α–position.
Heuristics for this improved version:
  1. The molecule must contain rings. We extract all ring atom indices.
  2. We build a graph among rings, connecting rings that share at least two atoms.
  3. We select the largest connected (fused) ring set and require that it is exactly 4 rings.
  4. For a classical steroid nucleus the four rings should be composed entirely of carbons,
     and their sizes must be three rings of 6 atoms and one ring of 5 atoms.
  5. Finally, we check that at least one [C@@H](O) substructure (an alpha–OH proxy)
     has its carbon atom residing in the fused nucleus.
If all these conditions are met, we classify the molecule as a 3α–hydroxy steroid.
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_3alpha_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 3α–hydroxy steroid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule qualifies as a 3α–hydroxy steroid, False otherwise.
        str: Explanation for the decision.
    """
    # Parse SMILES and assign stereochemistry.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    
    # Get ring information.
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    if not atom_rings:
        return False, "No rings found, not a steroid"
    
    # Build a graph of rings: each node is a ring (its index in atom_rings).
    # Connect two rings if they share at least two atoms (indicating fused rings).
    n_rings = len(atom_rings)
    adj = {i: set() for i in range(n_rings)}
    for i in range(n_rings):
        for j in range(i + 1, n_rings):
            if len(set(atom_rings[i]).intersection(atom_rings[j])) >= 2:
                adj[i].add(j)
                adj[j].add(i)
    
    # Identify connected components (fused systems) using depth-first search.
    seen = set()
    fused_components = []
    for i in range(n_rings):
        if i in seen:
            continue
        stack = [i]
        comp = set()
        while stack:
            node = stack.pop()
            if node in comp:
                continue
            comp.add(node)
            for neigh in adj[node]:
                if neigh not in comp:
                    stack.append(neigh)
        seen.update(comp)
        fused_components.append(comp)
    
    if not fused_components:
        return False, "No fused ring systems found, not a steroid nucleus"
    
    # Select largest fused ring system by number of rings.
    largest_component = max(fused_components, key=lambda comp: len(comp))
    
    # For classical steroids, we require exactly 4 fused rings.
    if len(largest_component) != 4:
        return False, f"Fused ring system has {len(largest_component)} rings (expected 4 for a steroid nucleus)"
    
    # Check that every ring in the chosen fused set is either 5- or 6-membered.
    ring_sizes = []
    for idx in largest_component:
        ring_size = len(atom_rings[idx])
        if ring_size not in (5, 6):
            return False, "Fused ring system contains rings that are not 5- or 6-membered"
        ring_sizes.append(ring_size)
    
    # For a classical steroid nucleus, we expect three 6-membered rings and one 5-membered ring.
    if sorted(ring_sizes) != sorted([5, 6, 6, 6]):
        return False, f"Fused ring system ring sizes {sorted(ring_sizes)} do not match classical steroid (5,6,6,6) pattern"
    
    # Build the set of all atom indices within the fused ring system.
    fused_atoms = set()
    for idx in largest_component:
        fused_atoms.update(atom_rings[idx])
    
    # Ensure that every atom in the fused nucleus is carbon (atomic number 6).
    for atom_idx in fused_atoms:
        if mol.GetAtomWithIdx(atom_idx).GetAtomicNum() != 6:
            return False, "Fused ring nucleus contains heteroatoms; expected an all–carbon steroid nucleus"
    
    # Look for an alpha–hydroxyl group by using SMARTS "[C@@H](O)".
    alpha_oh_smarts = "[C@@H](O)"
    alpha_oh_pattern = Chem.MolFromSmarts(alpha_oh_smarts)
    if alpha_oh_pattern is None:
        return False, "Invalid SMARTS for alpha hydroxyl group"
    
    matches = mol.GetSubstructMatches(alpha_oh_pattern)
    if not matches:
        return False, "No [C@@H](O) substructure found (alpha hydroxyl group not detected)"
    
    # Confirm that at least one match places the carbon (first atom in the match) in the fused nucleus.
    for match in matches:
        if match[0] in fused_atoms:
            return True, ("Molecule has a classical fused steroid nucleus (4 rings, 5+6+6+6, all carbons) "
                          "with a 3α-hydroxy substituent attached.")
    
    return False, "Alpha hydroxyl group not found on the fused steroid nucleus"

# Testing block:
if __name__ == "__main__":
    # Example: 5β-Cholestane-3α,7α,12α,23-tetrol
    test_smiles = "C1C[C@H](C[C@]2(C[C@H]([C@@]3([C@](C[C@@H]([C@]4([C@]3(CC[C@@]4([C@@](CC(CC(C)C)O)(C)[H])[H])[H])C)O)([C@@]12C)[H])[H])O)[H])O"
    result, reason = is_3alpha_hydroxy_steroid(test_smiles)
    print("Result:", result)
    print("Reason:", reason)