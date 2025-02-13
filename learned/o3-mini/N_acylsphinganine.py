"""
Classifies: CHEBI:31488 N-acylsphinganine
"""
"""
Classifies: N-acylsphinganine
A ceramide consisting of sphinganine in which one of the amino hydrogens is substituted by a fatty acyl group.
The improved version checks that the sphinganine backbone is present, and that the backbone nitrogen:
  - has exactly one hydrogen (it was originally NH2 and now becomes NH), 
  - is bonded to exactly two heavy atoms,
  - one neighbor is the backbone carbon and the other is a carbon that bears a carbonyl group.
"""

from rdkit import Chem

def is_N_acylsphinganine(smiles: str):
    """
    Determines if a molecule is an N-acylsphinganine based on its SMILES string.
    It checks for the presence of a sphinganine backbone fragment "N[C@@H](CO)[C@H](O)" 
    and then verifies that the backbone nitrogen is acylated.
    
    The criteria are:
      1. The molecule must contain the sphinganine core fragment.
      2. In that fragment the nitrogen must have exactly two heavy-atom neighbors.
         Furthermore, it should have exactly one hydrogen (consistent with N–acylation).
      3. The neighbor on the nitrogen that is not part of the sphinganine backbone must be a carbon atom.
      4. The bond from N to this acyl-candidate carbon must be single.
      5. This acyl candidate carbon must show at least one double bond to an oxygen (i.e. a carbonyl).
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule qualifies as an N-acylsphinganine, False otherwise.
        str: Explanation for the judgment.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS for the simple sphinganine backbone fragment.
    # This substructure approximates "N[C@@H](CO)[C@H](O)"
    sphinganine_pattern = Chem.MolFromSmarts("N[C@@H](CO)[C@H](O)")
    if sphinganine_pattern is None:
        return False, "Error in sphinganine SMARTS pattern"
    
    matches = mol.GetSubstructMatches(sphinganine_pattern)
    if not matches:
        return False, "Sphinganine backbone not found"
    
    # Iterate over each sphinganine match.
    for match in matches:
        # In the SMARTS the first atom is the nitrogen.
        n_idx = match[0]
        n_atom = mol.GetAtomWithIdx(n_idx)
        
        # Check that the nitrogen has exactly two heavy-atom neighbors.
        if n_atom.GetDegree() != 2:
            continue  # not the expected substitution pattern
        
        # In addition, require that the N carries exactly one hydrogen (N–acylation converts NH2 to NH)
        # Note: GetTotalNumHs counts both implicit and explicit Hs.
        if n_atom.GetTotalNumHs() != 1:
            continue
        
        # In our sphinganine SMARTS the second atom (index 1) is the alpha carbon of the backbone.
        backbone_neighbor_idx = match[1]
        
        # Identify the other neighbor, which should be the acyl group.
        acyl_candidate = None
        # Also check that the bond between N and the candidate is single.
        for bond in n_atom.GetBonds():
            nbr = bond.GetOtherAtom(n_atom)
            if nbr.GetIdx() == backbone_neighbor_idx:
                continue  # skip the backbone partner
            # Must be a single bond from N to acyl candidate.
            if bond.GetBondType() != Chem.BondType.SINGLE:
                continue
            acyl_candidate = nbr
            break
        
        if acyl_candidate is None:
            continue  # no valid acyl candidate found
        
        # The acyl candidate must be a carbon atom.
        if acyl_candidate.GetAtomicNum() != 6:
            continue
        
        # Check that the acyl candidate carbon is connected (via a double bond) to an oxygen (i.e. it bears a carbonyl).
        carbonyl_found = False
        for bond in acyl_candidate.GetBonds():
            # Skip the bond that goes back to the nitrogen.
            if bond.GetOtherAtom(acyl_candidate).GetIdx() == n_idx:
                continue
            # Look for a double bond to an oxygen atom.
            if bond.GetBondType() == Chem.BondType.DOUBLE:
                other = bond.GetOtherAtom(acyl_candidate)
                if other.GetAtomicNum() == 8:
                    carbonyl_found = True
                    break
        
        if not carbonyl_found:
            continue
        
        # (Optional) One might add additional rules here—for example, verifying that the acyl chain is long enough
        # or that the acyl candidate appears only once in the molecule—so as to avoid false positives from other
        # complicated structures.
        
        return True, "Molecule contains sphinganine backbone with proper N-acyl substitution"
    
    return False, "Sphinganine backbone found but nitrogen is not properly acylated (substitution pattern not matching)"

# For quick testing when running as a script.
if __name__ == "__main__":
    # Testing one of the true positive examples.
    test_smiles = "CCCCCCCCCCCCC[C@@H](O)[C@H](CO)NC(=O)CCCCCCCCC\\C=C/CCCCCCCC"
    result, reason = is_N_acylsphinganine(test_smiles)
    print("Test molecule:", result, reason)