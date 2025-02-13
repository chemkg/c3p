"""
Classifies: CHEBI:25106 macrolide
"""
"""
Classifies: A macrolide
Definition: A macrolide is a macrocyclic lactone with a ring of twelve or more atoms (derived from a polyketide).
Our approach is:
  1. Parse the molecule from its SMILES.
  2. Retrieve the ring information (SSSR rings).
  3. Identify “macrocyclic” rings with ≥12 atoms.
  4. Use an improved SMARTS pattern "[C;R](=O)[O;R]" so that both the carbonyl carbon and the ester oxygen are required to be in a ring.
  5. For each lactone candidate, check whether it is fully embedded in one of the macrocyclic rings.
  
Note: This scheme may still fail for highly fused or complex macrocycles, but it should improve over the previous attempt.
"""

from rdkit import Chem

def is_macrolide(smiles: str):
    """
    Determines if a molecule is a macrolide based on its SMILES string.
    A macrolide is defined as a macrocyclic lactone with a ring of 12 or more atoms.
    (The “derived from a polyketide” aspect is not explicitly tested.)

    The steps are:
      1. Parse the SMILES into a molecule.
      2. Identify all rings from the molecule (using SSSR information).
      3. From these, pick out rings with 12 or more atoms.
      4. Define a SMARTS for a lactone group that requires both the carbonyl carbon and the ester oxygen to be in a ring.
      5. Check if at least one match of this pattern is fully contained within a macrocyclic ring.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a macrolide, False otherwise.
        str: A text string providing the reason for the classification decision.
    """
    # Parse the molecule from SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Retrieve the SSSR ring information (as a list of tuples of atom indices)
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    
    # Identify macrocyclic rings (those with 12 or more atoms)
    macrocycle_rings = [set(ring) for ring in atom_rings if len(ring) >= 12]
    if not macrocycle_rings:
        return False, "No macrocyclic ring (12 or more atoms) found"

    # Define a SMARTS pattern for a lactone group.
    # By using "[C;R](=O)[O;R]", we require that both the carbonyl carbon (C;R) and the oxygen (O;R)
    # are themselves part of some ring. This helps focus on cyclic esters (lactones).
    lactone_pattern = Chem.MolFromSmarts("[C;R](=O)[O;R]")
    if lactone_pattern is None:
        return False, "Error in creating lactone SMARTS pattern"
    
    # Find all substructure matches for the lactone pattern.
    lactone_matches = mol.GetSubstructMatches(lactone_pattern)
    if not lactone_matches:
        return False, "No lactone (cyclic ester) group found"

    # For each lactone candidate, check if its atoms are contained in any of the macrocyclic rings.
    for match in lactone_matches:
        match_set = set(match)  # Typically contains (carbonyl carbon, oxygen)
        for ring in macrocycle_rings:
            if match_set.issubset(ring):
                return True, ("Found a macrolide: a macrocyclic ring (12+ atoms) contains "
                              "a lactone group ([C;R](=O)[O;R]).")
    
    # If lactone groups are found but none lie within a macrocycle
    return False, ("Found lactone group(s) but none are embedded in a macrocyclic (12 or more atoms) ring.")

# (Optional) Testing section:
if __name__ == "__main__":
    # Test examples: one expected to be a macrolide and one not.
    test_smiles_list = [
        # Paecilomycin K (should be a macrolide)
        "O=C1O[C@H](CC=C[C@@H](O)[C@@H]2O[C@H]([C@@H](C=3C1=C(O)C=C(OC)C3)O)CC2)C",
        # A small cyclic ester that is not a macrolide (ring size 4)
        "O=C1OC(C)C1"
    ]
    for smi in test_smiles_list:
        result, reason = is_macrolide(smi)
        print("SMILES:", smi)
        print("Result:", result)
        print("Reason:", reason)
        print("---")