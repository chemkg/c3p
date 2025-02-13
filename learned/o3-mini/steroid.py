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
  2. Secondary analysis: Find all rings of size 5 and 6, group fused rings (sharing >=2 atoms),
     and then look for a fused tetracyclic (4-ring) system that contains exactly one five-membered ring and three six-membered rings.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_steroid(smiles: str):
    """
    Determines if a molecule is a steroid based on its SMILES string.
    
    Uses two strategies:
      1. Primary substructure search using a SMARTS pattern for the canonical cyclopenta[a]phenanthrene nucleus.
      2. Secondary analysis: checking for a fused four-ring system (tetracycle) composed of exactly one 5-membered 
         ring and three 6-membered rings. This grouping is done by grouping rings (of size 5 or 6) into fused components 
         (two rings share 2+ atoms).
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a steroid, False otherwise.
        str: Explanation of the classification decision.
    """
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Primary check: Use a SMARTS pattern for the steroid nucleus.
    # This pattern represents the classic cyclopenta[a]phenanthrene core.
    steroid_pattern = Chem.MolFromSmarts("C1CCC2C3CCC4C1C3CCC2C4")
    if mol.HasSubstructMatch(steroid_pattern):
        return True, "Molecule contains the canonical steroid nucleus (cyclopenta[a]phenanthrene skeleton)."
    
    # Secondary check: Analyze the fused ring system.
    ring_info = mol.GetRingInfo()
    all_rings = ring_info.AtomRings()  # List of tuples (atom indices)
    
    # Filter rings: only consider rings of size 5 or 6
    candidate_rings = []
    for ring in all_rings:
        if len(ring) in (5, 6):
            candidate_rings.append(set(ring))
    
    if not candidate_rings:
        return False, "No 5- or 6-membered rings found in the molecule."
    
    # Build connectivity (fused rings) between these candidate rings:
    # Two rings are considered fused if they share at least 2 atoms.
    n = len(candidate_rings)
    adj_list = {i: [] for i in range(n)}
    for i in range(n):
        for j in range(i+1, n):
            # Compute the intersection size
            if len(candidate_rings[i].intersection(candidate_rings[j])) >= 2:
                adj_list[i].append(j)
                adj_list[j].append(i)
    
    # Find connected components among candidate rings using DFS.
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
    
    # Look for a fused component that is tetracyclic (4 rings) with exactly 1 five-membered and 3 six-membered rings.
    for comp in fused_components:
        if len(comp) == 4:
            count5 = 0
            count6 = 0
            for idx in comp:
                ring_size = len(candidate_rings[idx])
                if ring_size == 5:
                    count5 += 1
                elif ring_size == 6:
                    count6 += 1
            if count5 == 1 and count6 == 3:
                return True, ("Molecule contains a fused tetracyclic ring system (1 five-membered and 3 six-membered rings) "
                              "suggestive of a steroid nucleus (possibly modified).")
    
    return False, "Molecule does not appear to contain a steroid-like fused tetracyclic system."

# Example test cases (uncomment the following for interactive testing):
# test_smiles = [
#     "O=C1C[C@]2([C@@]([C@@]3(C([C@]4([C@@]([C@](CC4)([C@@H](CC/C(/C(C)C)=C/C)C)[H])(CC3)C)[H])=CC2)[H])(CC1)C",  # Avenastenone (steroid)
#     "S1S[C@]23N([C@@H]4[C@@H](O)[C@@H](O)CC([C@@H]4C2)=O)C([C@]15N([C@@H]6[C@@H](OC)[C@@H](O)CC([C@@H]6C5)=O)C3=O)=O",  # Epicoccin U (non-steroid)
# ]
# for s in test_smiles:
#     result, reason = is_steroid(s)
#     print("Result:", result)
#     print("Reason:", reason)
#     print()