"""
Classifies: CHEBI:83970 cardiac glycoside
"""
"""
Classifies: Cardiac glycoside (Steroid lactones containing sugar residues that act on the contractile force of the cardiac muscles)

Heuristic approach:
  1. Detect a steroid nucleus by finding a fused system of at least 4 rings.
  2. Detect a lactone (butenolide) ring – here using a simplified SMARTS.
  3. Detect at least one sugar residue (pyranose or furanose-like ring).
  
Note: The SMARTS patterns used and the ring-clustering algorithm are approximate.
"""

from rdkit import Chem

def is_cardiac_glycoside(smiles: str):
    """
    Determines if a molecule is a cardiac glycoside based on its SMILES string.
    Cardiac glycosides are defined as steroid lactones containing sugar residues.
    
    Heuristic criteria:
      - A steroid nucleus: a fused ring system with at least 4 rings that share bonds.
      - A lactone ring: a five‐membered unsaturated ring with a carbonyl and an oxygen (butenolide).
      - A sugar residue: detected using common pyranose or furanose ring SMARTS.
      
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule qualifies as a cardiac glycoside, False otherwise.
        str: Reason for the classification.
    """
    
    # Parse the input SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # --- Helper: check for fused steroid nucleus ---
    def has_steroid_nucleus(mol):
        """
        Determines if the molecule has a fused ring system with at least 4 rings.
        Approach:
          - Get all rings from the molecule.
          - Build a graph connecting rings that share at least 2 atoms (typical for fused rings).
          - If any connected component contains 4 or more rings, we count that as a steroid nucleus.
        """
        ring_info = mol.GetRingInfo()
        rings = ring_info.AtomRings()
        if len(rings) < 4:
            return False
        # Build a graph where each ring is a node.
        n = len(rings)
        adjacency = {i: set() for i in range(n)}
        for i in range(n):
            for j in range(i+1, n):
                # Two rings are “fused” if they share at least 2 atoms.
                if len(set(rings[i]).intersection(rings[j])) >= 2:
                    adjacency[i].add(j)
                    adjacency[j].add(i)
        # Find connected components in the ring graph.
        visited = set()
        def dfs(node, comp):
            comp.add(node)
            for neighbor in adjacency[node]:
                if neighbor not in comp:
                    dfs(neighbor, comp)
            return comp
        
        for i in range(n):
            if i not in visited:
                comp = dfs(i, set())
                if len(comp) >= 4:
                    return True
                visited.update(comp)
        return False
    
    if not has_steroid_nucleus(mol):
        return False, "No steroid nucleus (fused system of 4 or more rings) detected"
    
    # --- Check for lactone ring (butenolide) ---
    # Use a SMARTS that is meant to capture a five-membered lactone (furanone) moiety.
    # Note that many cardiac glycosides feature such a moiety.
    lactone_smarts = "C1=CC(=O)O1"  # simplified butenolide pattern
    lactone_pattern = Chem.MolFromSmarts(lactone_smarts)
    if not mol.HasSubstructMatch(lactone_pattern):
        return False, "No lactone (butenolide) ring detected"
    
    # --- Check for sugar residues ---
    # Look for signatures of pyranose or furanose rings.
    # The following SMARTS are simplified representations.
    pyranose_smarts = "OC1C(O)C(O)C(O)C(O)C1"   # typical six-membered sugar pattern (simplified)
    furanose_smarts = "OC1C(O)C(O)OC1"             # typical five-membered sugar pattern (simplified)
    pyranose_pattern = Chem.MolFromSmarts(pyranose_smarts)
    furanose_pattern = Chem.MolFromSmarts(furanose_smarts)
    
    sugar_found = False
    if mol.HasSubstructMatch(pyranose_pattern):
        sugar_found = True
    elif mol.HasSubstructMatch(furanose_pattern):
        sugar_found = True
        
    if not sugar_found:
        return False, "No sugar residue (glycosidic moiety) detected"
    
    # If all criteria are met, classify as cardiac glycoside.
    msg = ("Detected a fused steroid nucleus (≥4 fused rings), a lactone (butenolide) ring, "
           "and at least one sugar residue consistent with a cardiac glycoside structure.")
    return True, msg


# Example usage:
if __name__ == "__main__":
    # Test example with one of the provided SMILES (Erychroside).
    test_smiles = "O[C@@]12[C@]3([C@@]([C@@]4([C@](O)(CC3)C[C@@H](O[C@@H]5O[C@@H]([C@@H](O[C@@H]6OC[C@@H](O)[C@H](O)[C@H]6O)[C@@H](O)C5)C)CC4)C=O)(CC[C@@]1([C@H](CC2)C=7COC(=O)C7)C)[H])[H]"
    result, reason = is_cardiac_glycoside(test_smiles)
    print("Is cardiac glycoside?", result)
    print("Reason:", reason)