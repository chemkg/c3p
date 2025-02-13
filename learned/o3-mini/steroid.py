"""
Classifies: CHEBI:35341 steroid
"""
#!/usr/bin/env python
"""
Classifies: Steroid
Definition:
  Any of naturally occurring compounds and synthetic analogues, based on the cyclopenta[a]phenanthrene carbon skeleton,
  partially or completely hydrogenated; there are usually methyl groups at C-10 and C-13, and often an alkyl group at C-17.
  By extension, one or more bond scissions, ring expansions and/or ring contractions of the skeleton may have occurred.

Revised strategy:
  1. Primary SMARTS search for the canonical steroid nucleus.
  2. Secondary analysis: Find all 5- and 6-membered rings, then group fused rings (rings sharing >=2 atoms)
     and search for a fused tetracyclic component containing exactly 1 five-membered ring and 3 six-membered rings.
  3. Verify that the fused core has a total carbon count roughly matching a steroid nucleus (15–18 carbons).
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_steroid(smiles: str):
    """
    Determines if a molecule is a steroid based on its SMILES string.
    
    Two strategies are used:
      1. Primary: A direct substructure search using a SMARTS pattern for the classic cyclopenta[a]phenanthrene core.
      2. Secondary: Analyze the fused ring system for 5- and 6-membered rings.
         Look for a fused tetracyclic component (4 rings) having exactly one 5-membered and three 6-membered rings.
         Then, combine the atoms from these rings and check that the number of carbon atoms in the core is about 17.
         (We allow a range of 15–18 to account for minor modifications.)
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a steroid, False otherwise.
        str: Explanation of the classification decision.
    """
    # Parse the SMILES string using RDKit
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Primary check: SMARTS pattern for the canonical cyclopenta[a]phenanthrene skeleton
    # (This pattern represents four fused rings in the classical arrangement.)
    steroid_pattern = Chem.MolFromSmarts("C1CCC2C3CCC4C1C3CCC2C4")
    if mol.HasSubstructMatch(steroid_pattern):
        return True, "Molecule contains the canonical steroid nucleus (cyclopenta[a]phenanthrene skeleton)."
    
    # Secondary check: Analyze fused ring systems from rings of size 5 and 6
    ring_info = mol.GetRingInfo()
    all_rings = ring_info.AtomRings()  # List of tuples of atom indices for each ring
    
    # Only consider rings of size 5 or 6
    candidate_rings = []
    for ring in all_rings:
        if len(ring) in (5, 6):
            candidate_rings.append(set(ring))
    
    if not candidate_rings:
        return False, "No 5- or 6-membered rings found in the molecule."
    
    # Build connectivity between candidate rings: two rings are fused if they share ≥2 atoms.
    n = len(candidate_rings)
    adj_list = {i: [] for i in range(n)}
    for i in range(n):
        for j in range(i+1, n):
            if len(candidate_rings[i].intersection(candidate_rings[j])) >= 2:
                adj_list[i].append(j)
                adj_list[j].append(i)
    
    # Find connected components (fused ring groups) using DFS.
    visited = set()
    fused_components = []
    for i in range(n):
        if i not in visited:
            stack = [i]
            comp = set()
            while stack:
                node = stack.pop()
                if node in comp:
                    continue
                comp.add(node)
                visited.add(node)
                for neighbor in adj_list[node]:
                    if neighbor not in comp:
                        stack.append(neighbor)
            fused_components.append(comp)
    
    # Look for a fused component with 4 rings: expecting exactly 1 five-membered and 3 six-membered rings.
    for comp in fused_components:
        if len(comp) == 4:
            count5 = 0
            count6 = 0
            union_atoms = set()
            for idx in comp:
                ring_size = len(candidate_rings[idx])
                union_atoms |= candidate_rings[idx]
                if ring_size == 5:
                    count5 += 1
                elif ring_size == 6:
                    count6 += 1
            if count5 == 1 and count6 == 3:
                # Further check: count carbon atoms in this fused core.
                c_count = sum(1 for i in union_atoms if mol.GetAtomWithIdx(i).GetAtomicNum() == 6)
                if 15 <= c_count <= 18:
                    return True, ("Molecule contains a fused tetracyclic ring system (1 five-membered and 3 six-membered rings) "
                                  "with a core carbon count of {} (suggestive of a steroid nucleus).".format(c_count))
    
    return False, "Molecule does not appear to contain a steroid-like fused tetracyclic system with appropriate carbon composition."

# Example testing (uncomment for interactive testing):
# test_smiles = [
#     "O=C1C[C@]2([C@@]([C@@]3(C([C@]4([C@@]([C@](CC4)([C@@H](CC/C(/C(C)C)=C/C)C)[H])(CC3)C)[H])=CC2)[H])(CC1)C",  # Avenastenone (steroid)
#     "O=C1N2[C@H]([C@@]3(O[C@](C(N3[C@H]1CC(C)C)=O)(NC(=O)[C@]4(O)C=C5C6=C7C(NC=C7C[C@H]5N(C4)C)=CC=C6)C(C)C)O)CCC2",  # 8alpha-hydroxy-alpha-ergokryptine (false positive)
# ]
# for s in test_smiles:
#     result, reason = is_steroid(s)
#     print("Result:", result)
#     print("Reason:", reason)
#     print()