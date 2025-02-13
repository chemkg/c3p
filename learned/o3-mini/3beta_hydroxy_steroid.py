"""
Classifies: CHEBI:36836 3beta-hydroxy steroid
"""
#!/usr/bin/env python
"""
Classifies: CHEBI: 3β-hydroxy steroid
Defined as: A 3-hydroxy steroid in which the 3-hydroxy substituent is in the β–position.
Note:
  Because determining the exact position (the 3–position) and β–orientation from SMILES is nontrivial,
  this heuristic:
    1. Checks for a fused tetracyclic (4-ring) system via connectivity analysis of rings.
    2. Searches for an exocyclic hydroxyl group (–OH) attached to a chiral ring-carbon 
       whose CIP code is “R” (as typically found in 3β–hydroxy steroids such as cholesterol derivatives).
  This approach is heuristic and may fail for some cases.
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

    # Add explicit hydrogens so we can better detect -OH groups
    mol = Chem.AddHs(mol)
    
    # Ensure stereochemical information is assigned (this will compute CIP codes)
    AllChem.AssignStereochemistry(mol, cleanIt=True, force=True)

    # Get ring information from the molecule
    ring_info = mol.GetRingInfo()
    rings = [set(ring) for ring in ring_info.AtomRings()]
    if not rings:
        return False, "No rings found; steroid nucleus expected to have four fused rings"

    # Build connectivity between rings: two rings share an edge if they have at least 2 atoms in common.
    n_rings = len(rings)
    ring_adj = {i: set() for i in range(n_rings)}
    for i in range(n_rings):
        for j in range(i+1, n_rings):
            if len(rings[i].intersection(rings[j])) >= 2:
                ring_adj[i].add(j)
                ring_adj[j].add(i)

    # Find connected components among rings (fused blocks)
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
                for neighbor in ring_adj[node]:
                    if neighbor not in comp:
                        stack.append(neighbor)
            visited.update(comp)
            fused_components.append(comp)

    # We expect the steroid nucleus to have at least 4 fused rings.
    if not any(len(component) >= 4 for component in fused_components):
        max_rings = max((len(component) for component in fused_components), default=0)
        return False, f"Fused tetracyclic ring system not found; largest fused component contains {max_rings} ring(s)"

    # Now search for candidate hydroxyl groups.
    candidate_found = False
    candidate_msg = ""
    # Loop over oxygen atoms; check for -OH groups.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 8:
            continue  # not an oxygen        
        # Check if this oxygen has at least one hydrogen attached (looking for -OH)
        if atom.GetTotalNumHs() < 1:
            continue
        
        # Get all neighboring heavy atoms
        heavy_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() > 1]
        if len(heavy_neighbors) != 1:
            continue  # -OH should be attached to one heavy atom
        
        neighbor = heavy_neighbors[0]
        # Require that the oxygen is exocyclic (not part of a ring) but its attached carbon is in a ring.
        if atom.IsInRing():
            continue
        if neighbor.GetAtomicNum() != 6 or not neighbor.IsInRing():
            continue
        
        # Check whether the neighbor carbon is chiral (has a CIP code property)
        if not neighbor.HasProp('_CIPCode'):
            continue  # Cannot determine stereochemistry without CIP code
        cip = neighbor.GetProp('_CIPCode')
        
        # For 3β-hydroxy steroids (e.g. cholesterol) the chiral carbon at position 3 often has CIP code "R".
        if cip == "R":
            candidate_found = True
            candidate_msg = (f"Steroid nucleus detected (largest fused component ≥4 rings) and candidate "
                             f"3β-hydroxy group found on carbon {neighbor.GetIdx()} with CIP code {cip}.")
            break

    if not candidate_found:
        return False, ("Steroid nucleus may be present, but no candidate 3β-hydroxy group "
                       "(exocyclic -OH on a chiral ring-carbon with CIP code 'R') was found")
    
    return True, candidate_msg

# Example test cases for module testing:
if __name__ == "__main__":
    # Example steroid that should be classified as a positive case.
    test_smiles_positive = "CC(C)CCC[C@@H](C)[C@H]1CC[C@H]2C3=C(CC[C@]12C)[C@@]1(C)CC[C@H](O)[C@@](C)(CO)[C@@H]1CC3"
    result, reason = is_3beta_hydroxy_steroid(test_smiles_positive)
    print("Positive Test:")
    print("SMILES:", test_smiles_positive)
    print("Result:", result)
    print("Reason:", reason)
    print()

    # An example negative test: a non-steroid molecule.
    test_smiles_negative = "O=C1C2=CC=CC=C2C(=O)N1"
    result, reason = is_3beta_hydroxy_steroid(test_smiles_negative)
    print("Negative Test:")
    print("SMILES:", test_smiles_negative)
    print("Result:", result)
    print("Reason:", reason)