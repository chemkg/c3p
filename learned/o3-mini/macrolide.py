"""
Classifies: CHEBI:25106 macrolide
"""
"""
Classifies: A macrolide
Definition: A macrocyclic lactone with a ring of twelve or more members derived from a polyketide.
This classifier checks if a molecule has at least one ring of 12+ atoms that fully contains a lactone group 
(i.e. an ester functional group with both carbonyl carbon and oxygen in that ring).
"""

from rdkit import Chem

def is_macrolide(smiles: str):
    """
    Determines if a molecule is a macrolide based on its SMILES string. 
    A macrolide is defined as a macrocyclic lactone with a ring of 12 or more atoms. 
    (The "derived from a polyketide" aspect is not explicitly tested.)
    
    The approach is:
      1. Parse the molecule and retrieve ring information.
      2. Identify macrocyclic rings (12 or more atoms).
      3. Find ester substructure matches (which are candidates for lactone groups).
      4. For each ester match, check whether both the carbonyl carbon and the ester oxygen
         are contained together in one of the macrocycle rings.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule is a macrolide, False otherwise.
        str: Reason for classification.
    """
    # Parse the molecule from SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get all rings (each is a tuple of atom indices)
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    
    # Identify macrocyclic rings: rings with 12 or more atoms
    macrocycle_rings = [set(ring) for ring in rings if len(ring) >= 12]
    if not macrocycle_rings:
        return False, "No macrocyclic ring (12 or more atoms) found"
    
    # Define a SMARTS pattern for an ester group.
    # This pattern matches a carbonyl carbon ([CX3](=O)) attached to an oxygen ([OX2]).
    # It may match esters that are not lactones, so we will check ring membership later.
    lactone_candidate_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2]")
    if lactone_candidate_pattern is None:
        return False, "Error in creating ester SMARTS pattern"
    
    # Get all ester substructure matches (as tuples of atom indices, e.g. (carbon, oxygen) )
    matches = mol.GetSubstructMatches(lactone_candidate_pattern)
    if not matches:
        return False, "No lactone (cyclic ester) group found"
    
    # For each ester candidate, check if both the carbonyl carbon and the ester oxygen
    # are found together in any macrocyclic ring.
    for match in matches:
        match_set = set(match)  # Typically (carbonyl carbon, oxygen)
        # Loop over each macrocycle ring
        for ring in macrocycle_rings:
            if match_set.issubset(ring):
                return True, "Found a macrolide: macrocyclic ring (12+ atoms) containing a lactone group."
    
    return False, "Found lactone group(s) but none are embedded in a macrocyclic (12 or more atoms) ring."

# (Optional) Testing section:
if __name__ == "__main__":
    # List one or two test SMILES (feel free to extend with any candidate from the provided list)
    test_smiles_list = [
        # Paecilomycin K (should be macrolide): 
        "O=C1O[C@H](CC=C[C@@H](O)[C@@H]2O[C@H]([C@@H](C=3C1=C(O)C=C(OC)C3)O)CC2)C",
        # A structure that is not macrolide: a small cyclic ester
        "O=C1OC(C)C1"
    ]
    for smi in test_smiles_list:
        result, reason = is_macrolide(smi)
        print("SMILES:", smi)
        print("Result:", result)
        print("Reason:", reason)
        print("---")