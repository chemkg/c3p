"""
Classifies: CHEBI:25106 macrolide
"""
"""
Classifies: A macrolide
Definition: A macrolide is a macrocyclic lactone with a ring of twelve or more atoms (derived from a polyketide).
Our strategy has been improved by merging SSSR rings that share atoms into larger macrocyclic candidates.
We then check if a lactone group (defined via SMARTS "[C;R](=O)[O;R]") is fully embedded in one of those macrocycles.
Note: The “derived from a polyketide” aspect is not explicitly tested.
"""

from rdkit import Chem

def is_macrolide(smiles: str):
    """
    Determines if a molecule is a macrolide based on its SMILES string.
    A macrolide is defined as a macrocyclic lactone with a cycle of 12 or more atoms.
    
    Steps:
      1. Parse the SMILES into a molecule.
      2. Retrieve SSSR ring information.
      3. Merge rings that share atoms into larger cycles (to account for fused rings).
      4. Select macrocycle candidates that have at least 12 atoms.
      5. Define a lactone SMARTS pattern "[C;R](=O)[O;R]" (cyclic ester, 
         requiring both the carbonyl carbon and the ester oxygen to be in a ring).
      6. Check if at least one lactone group is fully embedded in one of the macrocycle candidates.

    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a macrolide, False otherwise.
        str: Explanation of the classification decision.
    """
    # Parse the molecule from the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Retrieve SSSR ring information: each ring is a tuple of atom indices.
    ring_info = mol.GetRingInfo()
    simple_rings = list(ring_info.AtomRings())
    if not simple_rings:
        return False, "No rings found in molecule"
    
    # Convert each ring into a set of atom indices
    ring_sets = [set(ring) for ring in simple_rings]
    
    # Merge rings that share at least one atom; use a union-find like approach.
    n = len(ring_sets)
    parent = list(range(n))
    
    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i
    
    def union(i, j):
        ri = find(i)
        rj = find(j)
        if ri != rj:
            parent[rj] = ri
            
    # For every pair of rings, if they share atoms, union their components.
    for i in range(n):
        for j in range(i+1, n):
            if ring_sets[i].intersection(ring_sets[j]):
                union(i, j)
    
    # Group rings by their representative, and form the union of rings in each component.
    comp_dict = {}
    for i in range(n):
        root = find(i)
        if root not in comp_dict:
            comp_dict[root] = set()
        comp_dict[root].update(ring_sets[i])
    
    # Filter the unioned cycles to macrocycle candidates with 12 or more atoms.
    macrocycle_candidates = [atoms for atoms in comp_dict.values() if len(atoms) >= 12]
    if not macrocycle_candidates:
        return False, "No macrocyclic ring (12 or more atoms) found"
    
    # Define a SMARTS pattern for a lactone group:
    # "[C;R](=O)[O;R]" means the carbonyl carbon and the ester oxygen must both be in a ring.
    lactone_pattern = Chem.MolFromSmarts("[C;R](=O)[O;R]")
    if lactone_pattern is None:
        return False, "Error creating lactone SMARTS pattern"
    
    # Find all lactone group matches (each is a tuple of atom indices, usually (carbon, oxygen))
    lactone_matches = mol.GetSubstructMatches(lactone_pattern)
    if not lactone_matches:
        return False, "No lactone (cyclic ester) group found"
    
    # For each lactone match, check if it is fully embedded in a macrocyclic candidate.
    for match in lactone_matches:
        match_set = set(match)
        for macro_atoms in macrocycle_candidates:
            if match_set.issubset(macro_atoms):
                return True, ("Macrolide confirmed: found a lactone group that is embedded "
                              "within a macrocyclic ring containing {} atoms.".format(len(macro_atoms)))
    
    # If lactone groups are found but not contained entirely in any macrocycle candidate:
    return False, "Found lactone group(s) but none are embedded in a macrocyclic (12+ atoms) ring."

# (Optional) Example testing section:
if __name__ == "__main__":
    test_cases = [
        # Paecilomycin K (expected to be a macrolide)
        ("O=C1O[C@H](CC=C[C@@H](O)[C@@H]2O[C@H]([C@@H](C=3C1=C(O)C=C(OC)C3)O)CC2)C", True),
        # A small cyclic ester (not a macrolide, ring size too small)
        ("O=C1OC(C)C1", False)
    ]
    
    for smi, expected in test_cases:
        result, reason = is_macrolide(smi)
        print("SMILES:", smi)
        print("Expected:", expected, "Got:", result)
        print("Reason:", reason)
        print("-----")