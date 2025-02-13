"""
Classifies: CHEBI:36836 3beta-hydroxy steroid
"""
#!/usr/bin/env python
"""
Classifies: CHEBI: 3β-hydroxy steroid
Defined as: A 3-hydroxy steroid in which the 3–hydroxy substituent is in the β–position.
Heuristic details:
  1. Detect a fused tetracyclic (at least 4-ring) system.
  2. Within that nucleus look for an exocyclic hydroxyl group (-OH) attached to a ring carbon.
  3. Use RDKit’s chiral center detection (FindMolChiralCenters) to check that the ring carbon
     carrying the hydroxyl has a CIP label of "R" (characteristic of many 3β–hydroxy steroids).
Note: This approach is heuristic and may not pick up every 3β-hydroxy steroid.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_3beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 3β-hydroxy steroid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule appears to be a 3β-hydroxy steroid.
        str: Explanation for the decision.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Add explicit hydrogens to help in the identification of -OH groups.
    mol = Chem.AddHs(mol)
    
    # Assign stereochemistry (this computes CIP labels)
    AllChem.AssignStereochemistry(mol, cleanIt=True, force=True)
    
    # Get ring information; we expect a steroid nucleus to have four fused rings.
    ring_info = mol.GetRingInfo()
    rings = [set(ring) for ring in ring_info.AtomRings()]
    if not rings:
        return False, "No rings found; steroid nucleus expected to have four fused rings"
    
    # Build connectivity between rings: rings are adjacent if they share at least 2 atoms.
    n_rings = len(rings)
    ring_adj = {i: set() for i in range(n_rings)}
    for i in range(n_rings):
        for j in range(i+1, n_rings):
            if len(rings[i].intersection(rings[j])) >= 2:
                ring_adj[i].add(j)
                ring_adj[j].add(i)
    
    # Find fused (connected) ring components.
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
    
    # Determine the largest fused ring set.
    largest_component = max(fused_components, key=lambda comp: len(comp))
    if len(largest_component) < 4:
        return False, f"Fused tetracyclic ring system not found; largest fused component contains {len(largest_component)} ring(s)"
    
    # Gather all atom indices that belong to any ring in the largest fused component.
    steroid_ring_atoms = set()
    for comp in largest_component:
        steroid_ring_atoms.update(rings[comp])
    
    # Get chiral centers with their CIP labels using RDKit.
    # This returns a list of tuples: (atom_index, 'R' or 'S') for atoms with defined chirality.
    stereo_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    chiral_dict = {idx: label for idx, label in stereo_centers}
    
    # Now search for candidate hydroxyl groups attached to the steroid nucleus.
    candidate_found = False
    candidate_msg = ""
    for atom in mol.GetAtoms():
        # Look for oxygen atoms.
        if atom.GetAtomicNum() != 8:
            continue
        # Check if the oxygen appears to be in an -OH group (has at least one explicit hydrogen).
        # Note: GetTotalNumHs() includes implicit hydrogens if they were added.
        if atom.GetTotalNumHs() < 1:
            continue
        
        # Ensure the oxygen is exocyclic: do not consider oxygens already in a ring.
        if atom.IsInRing():
            continue
        
        # Get neighboring heavy atoms; for an -OH group ideally this is exactly one.
        heavy_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() > 1]
        if len(heavy_neighbors) != 1:
            continue
        neighbor = heavy_neighbors[0]
        
        # The candidate oxygen should be attached to a carbon in a ring.
        if neighbor.GetAtomicNum() != 6 or not neighbor.IsInRing():
            continue
        
        # Further require that the attached carbon belongs to the steroid nucleus.
        if neighbor.GetIdx() not in steroid_ring_atoms:
            continue
        
        # At this point, look up the chiral center assignment of the neighbor.
        if neighbor.GetIdx() not in chiral_dict:
            continue  # No chiral information on this atom
        cip = chiral_dict[neighbor.GetIdx()]
        
        # In many 3β-hydroxy steroids, the substituent (the –OH) appears on a chiral carbon with CIP "R".
        if cip == "R":
            candidate_found = True
            candidate_msg = (f"Steroid nucleus detected (largest fused component has {len(largest_component)} rings) and candidate "
                             f"3β–hydroxy group found on carbon {neighbor.GetIdx()} with CIP code {cip}.")
            break

    if not candidate_found:
        return False, ("Steroid nucleus may be present, but no candidate 3β–hydroxy group "
                       "(exocyclic -OH on a chiral ring-carbon with CIP code 'R') was found")
    
    return True, candidate_msg

# Module testing
if __name__ == "__main__":
    # Example positive: one representative steroid expected to be a 3β-hydroxy steroid.
    test_smiles_positive = "CC(C)CCC[C@@H](C)[C@H]1CC[C@H]2C3=C(CC[C@]12C)[C@@]1(C)CC[C@H](O)[C@@](C)(CO)[C@@H]1CC3"
    result, reason = is_3beta_hydroxy_steroid(test_smiles_positive)
    print("Positive Test:")
    print("SMILES:", test_smiles_positive)
    print("Result:", result)
    print("Reason:", reason)
    print()
    
    # Example negative test: a non-steroid molecule.
    test_smiles_negative = "O=C1C2=CC=CC=C2C(=O)N1"
    result, reason = is_3beta_hydroxy_steroid(test_smiles_negative)
    print("Negative Test:")
    print("SMILES:", test_smiles_negative)
    print("Result:", result)
    print("Reason:", reason)