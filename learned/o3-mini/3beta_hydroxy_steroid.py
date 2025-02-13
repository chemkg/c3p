"""
Classifies: CHEBI:36836 3beta-hydroxy steroid
"""
#!/usr/bin/env python
"""
Classifies: 3β-hydroxy steroid
Defined as: A 3-hydroxy steroid in which the 3-hydroxy substituent is in the β–position.
Note:
  Because determining the exact position (the 3–position) and β–orientation
  from SMILES is nontrivial without full numbering and 3D analysis, this heuristic:
    1. Checks for a fused tetracyclic (4-ring) system via a connectivity analysis of rings.
    2. Searches for an exocyclic hydroxyl group (–OH) that is attached to a chiral ring-carbon 
       that carries a CIP code “S” (commonly found for a 3β–OH).
  This approach is heuristic and may fail for unusual cases.
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_3beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 3β-hydroxy steroid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a 3β-hydroxy steroid, else False.
        str: Explanation of the classification decision.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Add explicit hydrogens to help detect -OH groups
    mol = Chem.AddHs(mol)
    # Ensure stereochemical information is assigned
    AllChem.AssignStereochemistry(mol, cleanIt=True, force=True)
    
    # Obtain ring information (list of rings as tuples of atom indices)
    ri = mol.GetRingInfo()
    rings = [set(r) for r in ri.AtomRings()]
    if not rings:
        return False, "No rings found; steroid nucleus expected to have four fused rings"
    
    # Build a graph of rings: connect two rings if they share at least 2 atoms.
    # Each node is the index of a ring (from the list 'rings').
    n_rings = len(rings)
    adj = {i: set() for i in range(n_rings)}
    for i in range(n_rings):
        for j in range(i+1, n_rings):
            if len(rings[i].intersection(rings[j])) >= 2:
                adj[i].add(j)
                adj[j].add(i)
    
    # Now, find connected components among the rings
    visited = set()
    fused_components = []
    for i in range(n_rings):
        if i not in visited:
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
            visited.update(comp)
            fused_components.append(comp)
    
    # Check if any connected component contains at least 4 rings.
    found_tetracyclic = any(len(comp) >= 4 for comp in fused_components)
    if not found_tetracyclic:
        return False, ("Fused tetracyclic ring system not found; "
                       f"largest fused component contains {max((len(comp) for comp in fused_components), default=0)} ring(s)")
    
    # Now search for candidate hydroxyl groups.
    # We loop over oxygen atoms that are likely to be –OH groups.
    candidate_found = False
    candidate_msg = ""
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 8:
            continue
        # Check if this oxygen has at least one hydrogen
        # (GetTotalNumHs gives count of connected explicit+implicit hydrogens)
        if atom.GetTotalNumHs() < 1:
            continue
        # Check neighbors: for an -OH, there should be exactly one heavy-atom neighbor
        heavy_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() > 1]
        if len(heavy_neighbors) != 1:
            continue
        nbr = heavy_neighbors[0]
        # We want the –OH to be exocyclic: oxygen should not be in a ring, but its attached carbon should be.
        if atom.IsInRing():
            continue
        if nbr.GetAtomicNum() != 6 or not nbr.IsInRing():
            continue
        
        # Check if the neighbor carbon has any chiral tag
        if not nbr.HasProp('_CIPCode'):
            continue  # Without the CIP code we cannot be sure of the stereochemistry.
        cip = nbr.GetProp('_CIPCode')
        if cip == "S":
            candidate_found = True
            candidate_msg = (f"Steroid nucleus found (largest fused component ≥4 rings) and candidate "
                             f"3β-hydroxy group detected on carbon {nbr.GetIdx()} with CIP code {cip}.")
            break

    if not candidate_found:
        return False, ("Steroid nucleus may be present, but no candidate 3β-hydroxy group "
                       "(exocyclic –OH on a chiral ring-carbon with CIP code 'S') was found")
    
    return True, candidate_msg

# For module testing:
if __name__ == "__main__":
    # Some test cases (1 example that should be classified as positive and one negative)
    test_smiles_positive = "CC(C)CCC[C@@H](C)[C@H]1CC[C@H]2C3=C(CC[C@]12C)[C@@]1(C)CC[C@H](O)[C@@](C)(CO)[C@@H]1CC3"
    result, reason = is_3beta_hydroxy_steroid(test_smiles_positive)
    print("Positive Test:")
    print("SMILES:", test_smiles_positive)
    print("Result:", result)
    print("Reason:", reason)
    print()

    test_smiles_negative = "O=C1C2=CC=CC=C2C(=O)N1"  # for example, a phthalimide (non-steroid)
    result, reason = is_3beta_hydroxy_steroid(test_smiles_negative)
    print("Negative Test:")
    print("SMILES:", test_smiles_negative)
    print("Result:", result)
    print("Reason:", reason)