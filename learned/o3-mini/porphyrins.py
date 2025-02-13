"""
Classifies: CHEBI:26214 porphyrins
"""
"""
Classifies: Porphyrins 
(Natural pigments containing a fundamental skeleton of four pyrrole nuclei united through the 
alpha‐positions by four methine groups to form a macrocyclic structure.)

Our improved approach is:
1. Parse the SMILES with RDKit.
2. Identify candidate pyrrole rings: these are 5‐membered rings with exactly one nitrogen.
3. Check the overall ring systems for a macrocyclic ring (typically 15–18 atoms) that contains exactly 4 nitrogen atoms.
4. If there are at least 4 distinct pyrrole candidate rings and a macrocyclic ring meeting the above criteria,
   the molecule is classified as a porphyrin.

Note: This is a heuristic approach and may not catch every nuance.
"""

from rdkit import Chem

def is_porphyrins(smiles: str):
    """
    Determines if a molecule is a porphyrin based on its SMILES string.
    
    Our improved heuristic:
      • Identify candidate pyrrole rings as 5-membered rings having exactly one nitrogen atom.
      • Identify a macrocycle ring (size 15 to 18) that has exactly 4 nitrogen atoms.
          This is meant to capture the overall porphyrin backbone (4 pyrrole rings connected by methine bridges).
      • Only if both criteria are met does the molecule get classified as a porphyrin.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a porphyrin, False otherwise.
        str: Explanation of the classification decision.
    """
    # Parse the SMILES:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get ring systems as tuples of atom indices:
    ring_info = mol.GetRingInfo().AtomRings()
    
    # Collect candidate pyrrole rings: 5-membered rings with exactly one nitrogen
    pyrrole_rings = set()
    for ring in ring_info:
        if len(ring) == 5:
            n_count = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 7)
            if n_count == 1:
                # using frozenset ensures unique sets (order independent)
                pyrrole_rings.add(frozenset(ring))
    
    num_candidate_pyrroles = len(pyrrole_rings)
    
    if num_candidate_pyrroles < 4:
        return False, f"Found only {num_candidate_pyrroles} candidate pyrrole rings (5-membered rings with one N), need at least 4 for a porphyrin."

    # Now search for a macrocycle ring that spans the overall porphyrin core.
    # Under the assumption that the four pyrrole rings are connected by four methine bridges,
    # the overall macrocycle ring should be roughly 16 atoms.
    macrocycle_found = False
    macrocycle_size = None
    for ring in ring_info:
        # Allowing some slack (15 to 18 atoms) to account for minor variations.
        if 15 <= len(ring) <= 18:
            # Count nitrogen atoms in this ring:
            n_in_ring = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 7)
            if n_in_ring == 4:
                macrocycle_found = True
                macrocycle_size = len(ring)
                break
    if not macrocycle_found:
        return False, "No macrocyclic ring (size 15-18 with exactly 4 N) found; the four pyrrole candidates are not connected in a typical porphyrin pattern."
    
    # If both conditions are met, classify as porphyrin.
    return True, (f"Found {num_candidate_pyrroles} candidate pyrrole rings and a macrocycle "
                  f"ring of size {macrocycle_size} with 4 nitrogen atoms, consistent with a porphyrin macrocycle.")

# Example usage (the testing below is only demonstrative):
if __name__ == '__main__':
    # A test with a known porphyrin-like compound.
    # (This is only a sample SMILES and may not represent a real porphyrin derivative.)
    test_smiles = "C1=2N3C(C=C4N5=C(C=C6N7C=8C(=C9N(C=1)C)[Mg]0)C(=C6C)CCC(O)=O)C(=C4C)C=C" 
    result, reason = is_porphyrins(test_smiles)
    print("Is porphyrin?", result)
    print("Reason:", reason)