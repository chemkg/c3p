"""
Classifies: CHEBI:74716 withanolide
"""
"""
Classifies: Withanolide – any steroid lactone that is a C28 (or near) steroid with a modified side chain 
forming a lactone ring and its substituted derivatives.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_withanolide(smiles: str):
    """
    Determines if a molecule is a withanolide based on its SMILES string.
    Withanolides are defined as steroid lactones having a near C28 steroid nucleus with a modified 
    side chain that forms a lactone ring. Our heuristic requirements are:
       - Valid SMILES.
       - At least 27 carbon atoms.
       - A fused ring system (steroid nucleus) formed by at least four rings that largely
         conforms to the classic steroid pattern (one five membered and mostly six membered rings).
       - At least one lactone ring (cyclic ester) present.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule is classified as a withanolide, False otherwise.
        str: Reason for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Count carbon atoms.
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 27:
        return False, f"Only {c_count} carbon atoms detected, fewer than expected for a steroid nucleus"
    
    # Obtain ring information.
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()  # each ring is a tuple of atom indices
    if len(rings) < 4:
        return False, f"Only {len(rings)} rings found, expected at least 4 rings (steroid nucleus typically has 4 fused rings)"
    
    # Build a connectivity graph of rings.
    # Two rings are considered fused if they share at least two atoms.
    ring_neighbors = {i: set() for i in range(len(rings))}
    for i in range(len(rings)):
        for j in range(i+1, len(rings)):
            common = set(rings[i]).intersection(rings[j])
            if len(common) >= 2:
                ring_neighbors[i].add(j)
                ring_neighbors[j].add(i)
                
    # Find the largest connected component of rings (fused system) using DFS.
    visited = set()
    def dfs(i, component):
        component.add(i)
        for nb in ring_neighbors[i]:
            if nb not in component:
                dfs(nb, component)
        return component
    
    largest_component = set()
    for i in range(len(rings)):
        if i not in visited:
            comp = dfs(i, set())
            visited.update(comp)
            if len(comp) > len(largest_component):
                largest_component = comp
    if len(largest_component) < 4:
        return False, (f"Largest fused ring system has only {len(largest_component)} rings; "
                       "expected at least 4 rings forming a steroid nucleus")
    
    # Check that the fused ring system has a pattern roughly consistent with a steroid nucleus:
    # most common steroid nucleus (cyclopentanoperhydrophenanthrene) has one 5-membered ring and three 6-membered rings.
    five_count = 0
    six_count = 0
    for idx in largest_component:
        ring_size = len(rings[idx])
        if ring_size == 5:
            five_count += 1
        elif ring_size == 6:
            six_count += 1
    # Heuristic: at least one 5-membered ring and at least two 6-membered rings 
    # (allowing some variation in substituted derivatives).
    if five_count < 1 or six_count < 2:
        return False, (f"Fused ring system pattern not consistent with a steroid nucleus: "
                       f"found {five_count} five-membered and {six_count} six-membered rings")
    
    # Look for a lactone ring: a cyclic ester.
    # The SMARTS pattern looks for a carbonyl (C(=O)) bonded to an oxygen that is part of a ring.
    lactone_smarts = "[CX3](=O)[OX2r]"
    lactone_pattern = Chem.MolFromSmarts(lactone_smarts)
    if not mol.HasSubstructMatch(lactone_pattern):
        return False, "No lactone ring (cyclic ester) found"
    
    return True, ("Molecule has a near-steroid carbon count, a fused ring system with a pattern consistent with a steroid nucleus "
                  "and a lactone ring – consistent with a withanolide structure.")

# Example usage (uncomment to test):
# smiles_example = "C[C@@H]([C@H]1CC[C@H]2[C@@H]3C[C@H]4O[C@]44[C@@H](O)C=CC(=O)[C@]4(C)[C@H]3CC[C@]12C)[C@H]1CC(C)=C(CO)C(=O)O1"  # withaferin A-like structure
# result, reason = is_withanolide(smiles_example)
# print(result, reason)