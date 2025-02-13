"""
Classifies: CHEBI:26214 porphyrins
"""
"""
Classifies: Porphyrins
(Natural pigments containing a fundamental skeleton of four pyrrole nuclei united through 
the alpha-positions by four methine groups to form a macrocyclic structure.)

Our improved heuristic:
  1. Parse the SMILES.
  2. Identify “candidate pyrrole rings” (5‐membered rings having exactly one nitrogen atom).
  3. From these rings, extract the unique core nitrogen atom indices.
  4. Require that there are at least 4 candidate pyrrole rings and that exactly 4 unique nitrogen 
     atoms are present.
  5. Check whether one of the overall rings (macrocycles) has a size roughly in the range 15–20 
     atoms and contains all 4 of these candidate core nitrogen atoms.
  6. Only if both criteria are met, classify the molecule as a porphyrin.
  
Note:
  This is a heuristic and some real porphyrin derivatives (or non-porphyrins) might be missed 
  or misclassified. Further refinements (or substructure SMARTS queries) could be added.
"""

from rdkit import Chem

def is_porphyrins(smiles: str):
    """
    Determines if a molecule is a porphyrin based on its SMILES string.
    
    The heuristic:
      • Identify candidate pyrrole rings as 5-membered rings having exactly one nitrogen atom.
      • From these, extract the unique set of core nitrogen atoms.
      • Verify that there are at least 4 candidate pyrrole rings and exactly 4 unique core nitrogens.
      • Find a macrocyclic ring (of size 15 to 20 atoms) that contains all 4 nitrogen atoms.
      • Classify as porphyrin only if both conditions are met.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a porphyrin, False otherwise.
        str: Explanation of the classification decision.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get ring information from the molecule.
    ring_info = mol.GetRingInfo().AtomRings()
    
    # Find candidate pyrrole rings: 5-membered rings with exactly one nitrogen.
    candidate_pyrroles = []
    for ring in ring_info:
        if len(ring) == 5:
            n_count = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 7)
            if n_count == 1:
                candidate_pyrroles.append(set(ring))
                
    if len(candidate_pyrroles) < 4:
        return False, (f"Found only {len(candidate_pyrroles)} candidate pyrrole rings "
                       " (5-membered rings with one N), need at least 4 for a porphyrin.")
    
    # Extract the unique core nitrogen atoms from these candidate rings.
    core_nitrogens = set()
    for ring in candidate_pyrroles:
        for idx in ring:
            # Only count if the atom is nitrogen.
            if mol.GetAtomWithIdx(idx).GetAtomicNum() == 7:
                core_nitrogens.add(idx)
    if len(core_nitrogens) != 4:
        return False, (f"Found {len(core_nitrogens)} unique core nitrogen atoms from candidate "
                       "pyrrole rings, but exactly 4 are expected for a porphyrin.")
    
    # Now check for an overall macrocyclic ring that connects the four candidate pyrroles.
    # We allow a macrocycle ring in the size range 15 to 20 atoms.
    macrocycle_found = False
    macrocycle_size = None
    for ring in ring_info:
        if 15 <= len(ring) <= 20:
            # Check if all core nitrogen atoms are in this ring.
            if core_nitrogens.issubset(set(ring)):
                macrocycle_found = True
                macrocycle_size = len(ring)
                break
                
    if not macrocycle_found:
        return False, ("No macrocyclic ring (size 15-20 with all 4 core N) found; even though "
                       "candidate pyrrole rings exist, they are not connected by a typical porphyrin "
                       "macrocycle.")
    
    return True, (f"Found {len(candidate_pyrroles)} candidate pyrrole rings with core N atoms {core_nitrogens} "
                  f"and a macrocyclic ring of size {macrocycle_size} containing all 4 nitrogen atoms, "
                  "consistent with a porphyrin macrocycle.")

# Example usage (demonstration only):
if __name__ == '__main__':
    # Test with a sample SMILES (this one is for 7(1)-hydroxychlorophyllide a shown in the prompt).
    test_smiles = "C=12N3C(=CC4=NC(=CC=5N(C=6C(=C7N=C(C1)[C@H]([C@@H]7CCC(O)=O)C)[C@H](C(C6C5C)=O)C(=O)OC)[Mg]3)C(=C4CO)CC)C(=C2C)C=C"
    result, reason = is_porphyrins(test_smiles)
    print("Is porphyrin?", result)
    print("Reason:", reason)