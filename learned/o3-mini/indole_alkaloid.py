"""
Classifies: CHEBI:38958 indole alkaloid
"""
"""
Classifies: Indole alkaloid
Definition: An alkaloid containing an indole skeleton.
This implementation uses a SMARTS pattern to detect the indole core 
(c1ccc2c(c1)[nH]cc2) and then confirms that the matched atoms cover two fused rings:
a five-membered (pyrrole) ring and a six-membered (benzene) ring sharing exactly two atoms.
It also checks that the molecule overall contains nitrogen.
"""

from rdkit import Chem

def is_indole_alkaloid(smiles: str):
    """
    Determines if a given molecule is an indole alkaloid based on its SMILES string.
    The molecule qualifies if it contains an indole skeleton – which here is defined as 
    a fused bicyclic system with a five-membered aromatic (pyrrole, with exactly 1 nitrogen)
    ring fused (sharing exactly 2 atoms) with a six-membered aromatic (benzene) ring – 
    and if the molecule overall contains nitrogen.

    Args:
        smiles (str): The SMILES string of the molecule.

    Returns:
        bool: True if the molecule is classified as an indole alkaloid, False otherwise.
        str: A reason string explaining the classification decision.
    """
    # Parse the SMILES into a molecule object.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # First, check that the molecule contains at least one nitrogen atom.
    n_total = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    if n_total == 0:
        return False, "Molecule does not contain any nitrogen atoms; not an alkaloid"

    # Define a SMARTS pattern that describes the classical indole skeleton.
    # This pattern means: a benzene ring fused to a five-membered ring with an NH.
    indole_smarts = "c1ccc2c(c1)[nH]cc2"
    indole_query = Chem.MolFromSmarts(indole_smarts)
    if indole_query is None:
        return False, "Error in indole SMARTS pattern"
    
    # Get substructure matches for the indole pattern.
    matches = mol.GetSubstructMatches(indole_query)
    if not matches:
        return False, "No indole skeleton found"
    
    # Retrieve the ring information from the molecule:
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()  # list of tuples of atom indices for each ring
    
    valid_indole = False
    # For each match we found (each is a tuple of atom indices matching the SMARTS)
    # we try to “recover” the two rings that form an indole:
    for match in matches:
        match_set = set(match)
        # We expect the classical indole SMARTS to match 9 atoms in the full indole system.
        # However, depending on substitution the SMARTS match size may vary.
        # So we use the ring info: we look for one 5-membered and one 6-membered ring
        # whose atoms are a subset of the SMARTS match and share exactly two atoms.
        candidate_five = []
        candidate_six = []
        for ring in atom_rings:
            ring_set = set(ring)
            # Consider only rings fully inside the SMARTS match.
            if ring_set.issubset(match_set):
                if len(ring) == 5:
                    candidate_five.append(ring_set)
                elif len(ring) == 6:
                    candidate_six.append(ring_set)
        # Now check if any pair of one 5-membered ring and one 6-membered ring shares exactly 2 atoms.
        for ring5 in candidate_five:
            for ring6 in candidate_six:
                if len(ring5.intersection(ring6)) == 2:
                    valid_indole = True
                    break
            if valid_indole:
                break
        if valid_indole:
            break

    # If no valid fused ring system was found by the extra check, then report failure.
    if not valid_indole:
        return False, "No valid indole fused ring system found (5-membered ring fused with 6-membered ring)"
    
    return True, "Molecule contains an indole skeleton (5-membered aromatic ring fused with 6-membered aromatic ring) and nitrogen, consistent with an indole alkaloid"


# For testing (if run as a script), one may include test cases:
if __name__ == "__main__":
    # Example SMILES string for Ochropposinine
    test_smiles = "CC[C@]1(CN2CCC=3C4=CC(=C(C=C4NC3[C@@]2(C[C@@]1(CCO)[H])[H])OC)OC)[H]"
    result, reason = is_indole_alkaloid(test_smiles)
    print("Result:", result)
    print("Reason:", reason)