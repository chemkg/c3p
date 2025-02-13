"""
Classifies: CHEBI:46640 diketone
"""
"""
Classifies: A compound that contains exactly two ketone functionalities (diketone)

A diketone here is defined as a compound in which exactly two carbonyl (C=O) groups appear,
each with a carbon bound to two other carbons (thus excluding aldehydes and carbonyls bound to heteroatoms),
with no extra carbonyl functionality. In addition, we wish to avoid highly complex fused ring systems,
so we reject molecules that contain a large fused ring cluster (more than 3 fused rings).
"""

from rdkit import Chem
from rdkit.Chem import Descriptors

def is_diketone(smiles: str):
    """
    Determines if a molecule qualifies as a simple diketone based on its SMILES string.
    Our criteria require:
      - Exactly two ketone groups defined by the SMARTS: [#6][CX3](=O)[#6]
      - Exactly two carbonyl functionalities (matched with [CX3]=O)
      - Not containing a highly complex (fused) ring system (i.e. fused cluster of more than 3 rings)
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule qualifies as a (simple) diketone, False otherwise.
        str: Explanation of the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS:
    # Ketone group: a carbonyl where the carbon is bonded to two carbons.
    ketone_smarts = "[#6][CX3](=O)[#6]"
    ketone_pattern = Chem.MolFromSmarts(ketone_smarts)
    # General carbonyl pattern (matches any C(=O))
    carbonyl_smarts = "[CX3]=O"
    carbonyl_pattern = Chem.MolFromSmarts(carbonyl_smarts)

    ketone_matches = mol.GetSubstructMatches(ketone_pattern)
    carbonyl_matches = mol.GetSubstructMatches(carbonyl_pattern)

    n_ketones = len(ketone_matches)
    n_carbonyls = len(carbonyl_matches)

    # Check that there are exactly 2 ketone groups by our definition.
    if n_ketones != 2:
        return False, f"Compound contains {n_ketones} ketone group(s) (by strict SMARTS criteria), which does not equal 2."
    
    # Enforce that no extra carbonyl functionality is present.
    if n_carbonyls != 2:
        return False, f"Compound contains {n_carbonyls} carbonyl group(s); extra carbonyl functionality detected beyond the 2 ketone groups."
    
    # Heuristic to penalize overly complex (fused) ring systems:
    # We compute the fused ring clusters from the molecule.
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()  # each ring is a tuple of atom indices

    # Build a graph between rings if they share at least one atom.
    n_rings = len(atom_rings)
    if n_rings > 0:
        # Create an adjacency list for rings.
        ring_graph = {i: set() for i in range(n_rings)}
        for i in range(n_rings):
            for j in range(i+1, n_rings):
                # If the rings share any atom, consider them fused.
                if set(atom_rings[i]).intersection(atom_rings[j]):
                    ring_graph[i].add(j)
                    ring_graph[j].add(i)
        # Find connected components (clusters) of rings.
        visited = set()
        fused_cluster_sizes = []
        for i in range(n_rings):
            if i not in visited:
                stack = [i]
                component = set()
                while stack:
                    cur = stack.pop()
                    if cur in visited:
                        continue
                    visited.add(cur)
                    component.add(cur)
                    for neigh in ring_graph[cur]:
                        if neigh not in visited:
                            stack.append(neigh)
                fused_cluster_sizes.append(len(component))
        max_fused = max(fused_cluster_sizes) if fused_cluster_sizes else 0
    else:
        max_fused = 0

    # Reject compounds that contain a large fused ring system.
    if max_fused > 3:
        return False, f"Compound has a complex fused ring system (cluster of {max_fused} rings) not typical for a simple diketone."
    
    # (Optional) Additional heuristic: if there are rings present and the molecular weight is high,
    # it might be a polycyclic natural product. For example, if the molecule has 3+ rings and MW >350
    # we could also flag it. Uncomment the following if desired:
    #
    # mw = Descriptors.ExactMolWt(mol)
    # num_rings = ring_info.NumRings()
    # if num_rings >= 3 and mw > 350:
    #     return False, f"Compound is likely a complex polycyclic structure (MW = {mw:.1f} and {num_rings} rings) not typical for a simple diketone."

    return True, "Contains exactly 2 ketone (C=O) groups (with both substituents being carbon) and no extra carbonyl functionality, with no overly complex fused ring system."

# Example usage (for testing purposes only)
if __name__ == "__main__":
    test_cases = [
        # True positives
        ("O=C(CCCCCCCCCCCCCCC)CC(=O)CCCCCCCCC", "10,12-Heptacosanedione"),
        ("C1[C@H](C(C([C@]2([H])[C@]1(C)[C@@]3(C(C=C([C@H]([C@]3(CC2)[H])C)C=C)=O)[H])(C)C)=O)O", "(+)-phytocassane A"),
        ("[H][C-](C(C)=O)C(C)=O", "acetylacetonate"),
        ("CCC(C(CC)=O)=O", "3,4-Hexanedione"),
        # False positive example (expected to be rejected by fused ring heuristic)
        ("O=C1C=C2C(=COC(=C2)/C=C/C)[C@@H]([C@@]1(O)C)CC(=O)CCCCCCC", "Monascusazaphilone B"),
        # False negative example: extra carbonyl groups – should be rejected.
        ("C1(C(=O)C(=O)C(C(=O)C1=O)O)O", "3,6-dihydroxycyclohexane-1,2,4,5-tetrone")
    ]
    for smi, name in test_cases:
        result, reason = is_diketone(smi)
        print(f"SMILES: {smi}")
        print(f"Name: {name}")
        print(f"Result: {result}, Reason: {reason}")
        print("-" * 60)