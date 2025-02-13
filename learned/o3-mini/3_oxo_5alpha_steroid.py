"""
Classifies: CHEBI:13601 3-oxo-5alpha-steroid
"""
"""
Classifies: CHEBI: 3-oxo-5alpha-steroid
Definition: A 3-oxo steroid that has alpha configuration at position 5.
Heuristic (improved):
  (1) Confirm a ketone group (C=O) on a ring.
  (2) Look for a steroid nucleus defined as a fused system of 4 rings—
      exactly one five-membered ring and three six-membered rings—
      whose union comprises 17 carbon atoms.
  (3) Require that at least one chiral centre is explicitly specified
      in the atoms forming the nucleus (a crude surrogate for the 5α designation).
Note: This approach is heuristic and may still mis‐classify borderline cases.
"""

from rdkit import Chem
from itertools import combinations

def is_3_oxo_5alpha_steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-5alpha-steroid based on its SMILES string.
    Heuristic:
      - Must have a ketone group in a ring.
      - Must contain a fused steroid nucleus. Here we require a set of 4 rings
        (1 five-membered and 3 six-membered) that are connected (fused; sharing at least 2 atoms)
        whose union gives 17 carbon atoms.
      - At least one chiral center must be annotated in the nucleus.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a 3-oxo-5alpha-steroid, else False.
        str: Explanation for the classification.
    """
    # Parse SMILES into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # Step 1: Look for the ketone group.
    # We use SMARTS to detect a C(=O) where the carbon is in a ring.
    ketone_pattern = Chem.MolFromSmarts("[#6;R](=O)")
    if not mol.HasSubstructMatch(ketone_pattern):
        return False, "No ring ketone found (3-oxo requirement not met)."
    
    # Step 2: Identify candidate rings from the molecule.
    ring_info = mol.GetRingInfo()
    all_rings = ring_info.AtomRings()  # Each is a tuple of atom indices.
    # Consider only rings of size 5 or 6 (typical sizes in steroid nucleus).
    candidate_rings = [ring for ring in all_rings if len(ring) in (5,6)]
    if len(candidate_rings) < 4:
        return False, f"Only {len(candidate_rings)} rings of size 5 or 6 found (steroid nucleus not detected)."
    
    # Helper: Check whether a given set of rings (each a tuple of atom indices) forms a connected fused system.
    def rings_connected(rings_subset):
        # Build a graph where each ring is a node.
        n = len(rings_subset)
        adjacency = {i: set() for i in range(n)}
        for i in range(n):
            for j in range(i+1, n):
                # Define fusion as sharing at least 2 atoms.
                if len(set(rings_subset[i]).intersection(rings_subset[j])) >= 2:
                    adjacency[i].add(j)
                    adjacency[j].add(i)
        # Do a DFS over this ring-graph.
        visited = set()
        def dfs(i):
            visited.add(i)
            for j in adjacency[i]:
                if j not in visited:
                    dfs(j)
        dfs(0)
        return len(visited) == n

    # Step 3: Search for a set of 4 rings that looks like a steroid nucleus.
    nucleus_found = False
    nucleus_atoms = set()
    for subset in combinations(candidate_rings, 4):
        # Check that the 4 rings are fused (each ring connects with at least one other,
        # and overall they form a connected component).
        if not rings_connected(subset):
            continue
        
        # Check ring sizes: we expect one five-membered and three six-membered rings.
        ring_sizes = [len(r) for r in subset]
        if ring_sizes.count(5) != 1 or ring_sizes.count(6) != 3:
            continue
        
        # Get the union of atoms in these rings.
        union_atoms = set()
        for ring in subset:
            union_atoms.update(ring)
        
        # Count how many of these atoms are carbons.
        carbon_count = sum(1 for idx in union_atoms if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
        if carbon_count == 17:
            nucleus_found = True
            nucleus_atoms = union_atoms
            break
    if not nucleus_found:
        return False, "Steroid nucleus not detected: no fused set of 4 rings (1 five-membered and 3 six-membered) with 17 carbons found."
    
    # Step 4: Confirm that the ketone group is located within the steroid nucleus.
    ketone_matches = mol.GetSubstructMatches(ketone_pattern)
    if not any(set(match).issubset(nucleus_atoms) for match in ketone_matches):
        return False, "Ketone group not found within the steroid nucleus."
    
    # Step 5: Check for explicit chiral information in the nucleus.
    # We require that at least one atom in the nucleus is assigned a chiral tag.
    chiral_in_nucleus = False
    for idx in nucleus_atoms:
        atom = mol.GetAtomWithIdx(idx)
        if atom.GetChiralTag() != Chem.rdchem.ChiralType.CHI_UNSPECIFIED:
            chiral_in_nucleus = True
            break
    if not chiral_in_nucleus:
        return False, "No chiral centers detected within the nucleus (5α configuration not confirmed)."
    
    # Passed all tests.
    return True, "Molecule has a steroid nucleus (4 fused rings: one 5-membered and three 6-membered rings with 17 carbons), a ketone on the nucleus, and explicit chiral annotation (consistent with a 3-oxo-5α-steroid)."


# Example usage (for testing; remove or comment out in production):
if __name__ == "__main__":
    # Test with one example: 5alpha-Androstane-3,11,17-trione
    test_smiles = "C[C@]12CC(=O)[C@H]3[C@@H](CC[C@H]4CC(=O)CC[C@]34C)[C@@H]1CCC2=O"
    result, reason = is_3_oxo_5alpha_steroid(test_smiles)
    print(result, reason)