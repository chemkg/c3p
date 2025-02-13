"""
Classifies: CHEBI:31488 N-acylsphinganine
"""
"""
Classifies: N-acylsphinganine
A ceramide consisting of sphinganine in which one of the amino hydrogens is substituted by a fatty acyl group.
"""

from rdkit import Chem

def is_N_acylsphinganine(smiles: str):
    """
    Determines if a molecule is an N-acylsphinganine based on its SMILES string.
    The molecule is expected to have a sphinganine backbone (for example, a fragment like 
    "N[C@@H](CO)[C@H](O)") and the nitrogen of that backbone must be acylated, that is, 
    bound to a carbonyl (C(=O)) from a fatty acid.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is classified as an N-acylsphinganine, False otherwise
        str: Reason for classification
    """
    
    # Parse the input SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern for a sphinganine backbone fragment.
    # We are looking for a nitrogen atom that is directly attached to a chiral carbon with a CH2OH group,
    # which in turn is connected to another chiral carbon carrying an OH group.
    # This pattern approximates the sphinganine core "N[C@@H](CO)[C@H](O)".
    sphinganine_pattern = Chem.MolFromSmarts("N[C@@H](CO)[C@H](O)")
    if not sphinganine_pattern:
        return False, "Error in sphinganine SMARTS pattern"
    
    matches = mol.GetSubstructMatches(sphinganine_pattern)
    if not matches:
        return False, "Sphinganine backbone not found"
    
    # For each match, check if the nitrogen atom (first atom in the match) is acylated.
    # That means that besides being connected to the sphinganine carbon, it should also connect 
    # to a carbon that is part of a carbonyl (C=O) group.
    for match in matches:
        n_idx = match[0]  # index of the nitrogen in the sphinganine backbone
        n_atom = mol.GetAtomWithIdx(n_idx)
        # Get neighbors of the nitrogen; one neighbor is the sphinganine carbon (already in the match)
        # Additional neighbor(s) would come from acyl substitution.
        neighbors = n_atom.GetNeighbors()
        acyl_found = False
        for neighbor in neighbors:
            # Skip the neighbor involved in the sphinganine backbone (it is an sp3 carbon with OH)
            # Check if this neighboring atom is a carbonyl carbon.
            if neighbor.GetAtomicNum() != 6:
                continue  # not carbon
            # For each bond connected to the neighbor, check if one is a double bond to oxygen.
            for bond in neighbor.GetBonds():
                # Identify the other atom in the bond.
                other_atom = bond.GetOtherAtom(neighbor)
                if other_atom.GetAtomicNum() == 8 and bond.GetBondTypeAsDouble() == 2:
                    # We found a carbonyl group attached to the neighbor.
                    acyl_found = True
                    break
            if acyl_found:
                break
        if acyl_found:
            return True, "Molecule contains sphinganine backbone with N-acyl substitution"
    
    return False, "Sphinganine backbone found but its nitrogen is not acylated (no carbonyl connected)"

# The following if __name__ == "__main__": block can be used for quick testing.
if __name__ == "__main__":
    # Example: N-(11Z)-icosenoylsphinganine
    test_smiles = "CCCCCCCCCCCCC[C@@H](O)[C@H](CO)NC(=O)CCCCCCCCC\\C=C/CCCCCCCC"
    result, reason = is_N_acylsphinganine(test_smiles)
    print("Test molecule:", result, reason)