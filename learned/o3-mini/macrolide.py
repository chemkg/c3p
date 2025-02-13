"""
Classifies: CHEBI:25106 macrolide
"""
"""
Classifies: A macrolide
Definition: A macrocyclic lactone with a ring of twelve or more members derived from a polyketide.
This classifier checks if a molecule contains a ring with 12+ atoms and a lactone (cyclic ester)
group within one of those rings.
"""

from rdkit import Chem

def is_macrolide(smiles: str):
    """
    Determines if a molecule is a macrolide based on its SMILES string.
    A macrolide is defined as a macrocyclic lactone with a ring of 12 or more atoms.
    (Note: The "derived from a polyketide" aspect is not explicitly tested.)

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a macrolide, False otherwise.
        str: Reason for classification.
    """
    # Parse the SMILES string to get a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()  # list of tuples of atom indices
    
    # Identify macrocyclic rings: 12 or more atoms
    macrocycle_rings = [ring for ring in rings if len(ring) >= 12]
    if not macrocycle_rings:
        return False, "No macrocyclic ring (12 or more atoms) found"

    # Define a SMARTS pattern for a lactone (cyclic ester).
    # This pattern looks for a carbonyl carbon in a ring bonded to an oxygen,
    # which in turn is bonded to a carbon in a ring.
    lactone_pattern = Chem.MolFromSmarts("[C;R](=O)O[C;R]")
    if lactone_pattern is None:
        return False, "Error in creating lactone SMARTS pattern"
    
    # Get all matches for the lactone substructure
    lactone_matches = mol.GetSubstructMatches(lactone_pattern)
    if not lactone_matches:
        return False, "No lactone (cyclic ester) group found"

    # Check whether at least one lactone match is fully contained in one macrocyclic ring.
    for match in lactone_matches:
        match_set = set(match)
        for ring in macrocycle_rings:
            if match_set.issubset(set(ring)):
                return True, f"Found a macrolide: macrocyclic ring with 12+ atoms and embedded lactone group."
    
    return False, "Found lactone group(s) but none are embedded in a macrocyclic (12 or more atoms) ring."

# Below are some tests that you might run (disabled by default):
if __name__ == "__main__":
    # An example SMILES from the provided list. (You can replace with any candidate)
    test_smiles = "O=C1OC2C(C(OC3(C2)OC(C(C)CC3)CC(O)CC)CC(=O)C(C=CCCCC(C(CC(C(C(C(C=C1)C)O)C)O)O)(O)C)CC(C)C)C"
    is_macro, reason = is_macrolide(test_smiles)
    print(is_macro, reason)