"""
Classifies: CHEBI:31488 N-acylsphinganine
"""
"""
Classifies: N-acylsphinganine
A ceramide consisting of sphinganine in which one of the amino hydrogens is substituted by a fatty acyl group.
Improved version: The code now checks that the nitrogen in the sphinganine backbone has exactly two heavy-atom neighbors—
one coming from the sphinganine chain and the other from the acyl group (which must include a carbonyl).
"""

from rdkit import Chem

def is_N_acylsphinganine(smiles: str):
    """
    Determines if a molecule is an N-acylsphinganine based on its SMILES string.
    The molecule is expected to contain a sphinganine backbone fragment such as
    "N[C@@H](CO)[C@H](O)" and the nitrogen of that backbone must be acylated.
    In an N-acylsphinganine the nitrogen should be substituted by exactly one acyl group,
    meaning that its degree (number of heavy-atom neighbors) is 2:
    one attached to the sphinganine chain and the other attached to a carbonyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is classified as an N-acylsphinganine, False otherwise
        str: Reason for classification
    """
    
    # Parse the input SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern to detect the sphinganine core fragment.
    # This pattern approximates the sphinganine backbone "N[C@@H](CO)[C@H](O)".
    sphinganine_pattern = Chem.MolFromSmarts("N[C@@H](CO)[C@H](O)")
    if sphinganine_pattern is None:
        return False, "Error in sphinganine SMARTS pattern"
    
    matches = mol.GetSubstructMatches(sphinganine_pattern)
    if not matches:
        return False, "Sphinganine backbone not found"
    
    # Iterate over each sphinganine match.
    # In the matched pattern, the first atom (index 0) is the nitrogen.
    # In a proper N-acylsphinganine the N should have exactly 2 heavy-atom neighbors:
    # one in the sphinganine backbone and one from the acyl (fatty acid) substitution.
    for match in matches:
        n_idx = match[0]
        n_atom = mol.GetAtomWithIdx(n_idx)
        
        # Check that nitrogen has exactly 2 neighbors.
        if n_atom.GetDegree() != 2:
            continue  # extra substitutions are present, so skip this match
        
        # Identify the neighbor that is not part of the backbone.
        # In our sphinganine SMARTS the backbone carbon (the alpha-carbon) is match[1].
        backbone_neighbor_idx = match[1]
        
        # Get the two neighbors.
        neighbors = n_atom.GetNeighbors()
        acyl_candidate = None
        for neighbor in neighbors:
            if neighbor.GetIdx() != backbone_neighbor_idx:
                acyl_candidate = neighbor
                break
        
        if acyl_candidate is None:
            continue  # no candidate found
        
        # The acyl candidate should be a carbon atom.
        if acyl_candidate.GetAtomicNum() != 6:
            continue
        
        # Check that the acyl candidate has at least one double bond to oxygen (i.e. a carbonyl group).
        carbonyl_found = False
        for bond in acyl_candidate.GetBonds():
            # bond between candidate and connected atom
            other = bond.GetOtherAtom(acyl_candidate)
            if other.GetAtomicNum() == 8 and bond.GetBondType() == Chem.BondType.DOUBLE:
                carbonyl_found = True
                break
        if not carbonyl_found:
            continue
        
        # If we got here, we found a sphinganine backbone whose nitrogen is substituted by a carbonyl-bearing carbon.
        return True, "Molecule contains sphinganine backbone with proper N-acyl substitution"
    
    return False, "Sphinganine backbone found but its nitrogen is not properly acylated (wrong substitution pattern)"

# Quick testing when running as a script.
if __name__ == "__main__":
    # Example test: N-(11Z)-icosenoylsphinganine
    test_smiles = "CCCCCCCCCCCCC[C@@H](O)[C@H](CO)NC(=O)CCCCCCCCC\\C=C/CCCCCCCC"
    result, reason = is_N_acylsphinganine(test_smiles)
    print("Test molecule:", result, reason)