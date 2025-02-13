"""
Classifies: CHEBI:74716 withanolide
"""
"""
Classifies: Withanolide – any steroid lactone that is a C27/C28 steroid (or substituted derivative)
with a modified side chain forming a lactone ring.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_withanolide(smiles: str):
    """
    Determines if a molecule is a withanolide based on its SMILES string.
    A withanolide is defined as a steroid lactone having a C28 (or near C28) steroid nucleus
    with a modified side chain that forms a lactone ring. Our heuristic requirements are:
      - Valid SMILES.
      - At least 27 carbon atoms (to allow for close derivatives).
      - A fused ring system (steroid nucleus) formed by at least four rings.
        (Here rings are considered fused if they share at least two atoms.)
      - At least one lactone ring (cyclic ester) present.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule is classified as withanolide, False otherwise.
        str: Reason for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Count carbon atoms.
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 27:
        return False, f"Only {c_count} carbon atoms detected, fewer than expected for a C27/C28 steroid nucleus"
    
    # Obtain ring information.
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()  # each ring is a tuple of atom indices
    if len(rings) < 4:
        return False, f"Only {len(rings)} rings found, expected at least 4 rings (a fused steroid nucleus typically has 4)"
    
    # Build a connectivity graph of rings.
    # Two rings are connected (fused) if they share at least two atoms.
    ring_neighbors = {i: set() for i in range(len(rings))}
    for i in range(len(rings)):
        for j in range(i+1, len(rings)):
            # find common atoms
            common = set(rings[i]).intersection(rings[j])
            if len(common) >= 2:
                ring_neighbors[i].add(j)
                ring_neighbors[j].add(i)
    
    # Find the size of the largest connected component (using a simple DFS)
    visited = set()
    def dfs(i, component):
        component.add(i)
        for nb in ring_neighbors[i]:
            if nb not in component:
                dfs(nb, component)
        return component
    
    largest_component_size = 0
    for i in range(len(rings)):
        if i not in visited:
            comp = dfs(i, set())
            visited.update(comp)
            if len(comp) > largest_component_size:
                largest_component_size = len(comp)
    
    if largest_component_size < 4:
        return False, (f"Largest fused ring system has only {largest_component_size} rings; "
                       "expected at least 4 rings forming a steroid nucleus")
    
    # Look for lactone ring: a cyclic ester.
    # SMARTS below looks for a carbonyl (C(=O)) directly bonded to an oxygen (O) that is part of a ring.
    lactone_smarts = "[CX3](=O)[OX2r]"
    lactone_pattern = Chem.MolFromSmarts(lactone_smarts)
    if not mol.HasSubstructMatch(lactone_pattern):
        return False, "No lactone ring (cyclic ester) found"
    
    # If all tests passed, classify as withanolide.
    return True, ("Molecule has a near-steroid carbon count, a fused ring system with at least 4 rings "
                  "consistent with a steroid nucleus, and a lactone ring – consistent with a withanolide structure.")

# Example usage (uncomment to test):
# smiles_example = "C[C@@H]([C@H]1CC[C@H]2[C@@H]3C[C@H]4O[C@]44[C@@H](O)C=CC(=O)[C@]4(C)[C@H]3CC[C@]12C)[C@H]1CC(C)=C(CO)C(=O)O1"  # withaferin A-like structure
# result, reason = is_withanolide(smiles_example)
# print(result, reason)

# Note:
# This heuristic may still misclassify some edge cases; further refinement (e.g. using more detailed steroid nucleus SMARTS)
# could improve performance.