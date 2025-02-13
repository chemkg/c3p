"""
Classifies: CHEBI:13601 3-oxo-5alpha-steroid
"""
"""
Classifies: CHEBI: 3-oxo-5alpha-steroid
Definition: A 3-oxo steroid that has alpha configuration at position 5.
Heuristic (improved):
  (1) The molecule must contain at least one ketone group (C=O) on a ring.
  (2) The molecule must contain a steroid nucleus: a fused system of 4 rings (exactly one five-membered ring and three six-membered rings) composed solely of carbon atoms and whose union has exactly 17 carbons.
  (3) At least one atom in the nucleus must be marked with a chiral tag.
Note: This heuristic approach tries to limit combinatorial searches by first filtering candidate rings and then working only on connected ring systems.
"""

from rdkit import Chem
from itertools import combinations

def is_3_oxo_5alpha_steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-5alpha-steroid based on its SMILES string.
    Heuristics:
      - The molecule must have a ketone group (C=O) on a ring.
      - The molecule must contain a fused steroid nucleus: a set of 4 rings, exactly one five-membered
        and three six-membered rings. All atoms of these rings should be carbons and their union comprises exactly 17 carbons.
      - At least one atom in the nucleus must have explicit chiral annotation.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a 3-oxo-5alpha-steroid, else False.
        str: Explanation for the classification.
    """
    
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # Step 1: Look for a ketone group in a ring using SMARTS.
    # This SMARTS checks for a carbon in a ring double-bonded to oxygen.
    ketone_pattern = Chem.MolFromSmarts("[#6;R](=O)")
    ketone_matches = mol.GetSubstructMatches(ketone_pattern)
    if not ketone_matches:
        return False, "No ring ketone found (3-oxo requirement not met)."
    
    # Step 2: Obtain all rings in the molecule.
    ring_info = mol.GetRingInfo()
    all_rings = ring_info.AtomRings()  # Each ring is a tuple of atom indices.
    
    # Filter candidate rings: only rings of size 5 or 6 and composed entirely of carbons.
    candidate_rings = []
    for ring in all_rings:
        if len(ring) not in (5, 6):
            continue
        if all(mol.GetAtomWithIdx(i).GetAtomicNum() == 6 for i in ring):
            candidate_rings.append(ring)
    
    if len(candidate_rings) < 4:
        return False, f"Only {len(candidate_rings)} candidate rings (size 5 or 6 with carbon atoms only) found."
    
    # Step 3: Build a connectivity graph on candidate rings.
    # Two rings are "fused" if they share at least 2 atoms.
    n_rings = len(candidate_rings)
    ring_adj = {i: set() for i in range(n_rings)}
    for i in range(n_rings):
        for j in range(i+1, n_rings):
            if len(set(candidate_rings[i]).intersection(candidate_rings[j])) >= 2:
                ring_adj[i].add(j)
                ring_adj[j].add(i)
    
    # Find connected components (fused systems) of rings.
    components = []
    visited = set()
    for i in range(n_rings):
        if i in visited:
            continue
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
        visited |= comp
        components.append(comp)
    
    # Step 4: In each component, look for a subset of 4 rings that matches the steroid nucleus.
    nucleus_found = False
    nucleus_atoms = set()
    for comp in components:
        # Only consider components that have at least 4 rings.
        if len(comp) < 4:
            continue
        comp_rings = [candidate_rings[i] for i in comp]
        # To reduce combinations, if many rings exist in the component then check all 4-ring combinations within it.
        for subset in combinations(comp_rings, 4):
            # Check connectivity within this subset
            # (We require that the subset itself forms a connected graph via the earlier fusion criterion.)
            # We rebuild a small graph for just this subset.
            sub_adj = {k: set() for k in range(4)}
            for a in range(4):
                for b in range(a+1, 4):
                    if len(set(subset[a]).intersection(subset[b])) >= 2:
                        sub_adj[a].add(b)
                        sub_adj[b].add(a)
            connected = False
            visited_subset = set()
            def dfs(n):
                visited_subset.add(n)
                for m in sub_adj[n]:
                    if m not in visited_subset:
                        dfs(m)
            dfs(0)
            if len(visited_subset) != 4:
                continue
            
            # Check ring sizes: exactly one ring of size 5 and three rings of size 6.
            sizes = [len(ring) for ring in subset]
            if sizes.count(5) != 1 or sizes.count(6) != 3:
                continue
            
            # Determine the union of atoms in these 4 rings.
            union_atoms = set()
            for ring in subset:
                union_atoms.update(ring)
            # Count carbons in the nucleus (should be exactly 17).
            carbon_count = sum(1 for idx in union_atoms if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
            if carbon_count != 17:
                continue
            
            # Check that at least one ketone group is located within this nucleus.
            ketone_in_nucleus = False
            for match in ketone_matches:
                if set(match).issubset(union_atoms):
                    ketone_in_nucleus = True
                    break
            if not ketone_in_nucleus:
                continue
            
            # Check for explicit chiral tags in at least one atom in the nucleus.
            chiral_found = False
            for idx in union_atoms:
                atom = mol.GetAtomWithIdx(idx)
                if atom.GetChiralTag() != Chem.rdchem.ChiralType.CHI_UNSPECIFIED:
                    chiral_found = True
                    break
            if not chiral_found:
                continue
            
            # We have found a suitable steroid nucleus.
            nucleus_found = True
            nucleus_atoms = union_atoms
            break
        if nucleus_found:
            break

    if not nucleus_found:
        return False, "Steroid nucleus not detected: no fused set of 4 rings (1 five-membered and 3 six-membered rings with 17 carbons) having a ketone and chiral center found."
    
    return True, "Molecule has a steroid nucleus (4 fused rings: one 5-membered and three 6-membered rings with 17 carbons), a ketone on the nucleus, and explicit chiral annotation (consistent with a 3-oxo-5α-steroid)."

# Example usage (for testing; remove or comment out in production):
if __name__ == "__main__":
    # Test with one example: 5alpha-Androstane-3,11,17-trione
    test_smiles = "C[C@]12CC(=O)[C@H]3[C@@H](CC[C@H]4CC(=O)CC[C@]34C)[C@@H]1CCC2=O"
    result, reason = is_3_oxo_5alpha_steroid(test_smiles)
    print(result, reason)