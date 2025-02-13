"""
Classifies: CHEBI:31488 N-acylsphinganine
"""
"""
Classifies: N-acylsphinganine
A ceramide consisting of sphinganine in which one of the amino hydrogens is substituted by a fatty acyl group.
The classifier looks for the sphinganine backbone fragment and then requires that the backbone nitrogen:
  - has exactly two heavy-atom neighbors,
  - carries exactly one hydrogen (i.e. from NH2 to NH),
  - is bonded (via a single bond) to a carbon which in turn:
       * is a carbon atom (atomic number 6),
       * is sp2-hybridized,
       * bears exactly one double bond to an oxygen (i.e. a carbonyl).
This additional check helps avoid false positives from more complex ceramides.
"""
from rdkit import Chem

def is_N_acylsphinganine(smiles: str):
    """
    Determines if a molecule is an N-acylsphinganine based on its SMILES string.
    It first detects the sphinganine backbone fragment "N[C@@H](CO)[C@H](O)".
    Then, for each match it verifies that the backbone nitrogen is acylated:
      - The nitrogen has exactly two heavy-atom neighbors (one being the backbone carbon).
      - It carries exactly one hydrogen (consistent with a NH in an amide).
      - Its non-backbone neighbor is a carbon connected via a single bond.
      - That acyl candidate carbon is sp2 hybridized.
      - The acyl candidate carbon is double-bonded to exactly one oxygen (i.e. a carbonyl).
      
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule qualifies as an N-acylsphinganine, False otherwise.
        str: Explanation for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern to approximate the sphinganine backbone: 
    # "N[C@@H](CO)[C@H](O)"
    sphinganine_pattern = Chem.MolFromSmarts("N[C@@H](CO)[C@H](O)")
    if sphinganine_pattern is None:
        return False, "Error in sphinganine SMARTS pattern"
    
    matches = mol.GetSubstructMatches(sphinganine_pattern)
    if not matches:
        return False, "Sphinganine backbone not found"
    
    # Iterate over every sphinganine match.
    for match in matches:
        # The first atom in the pattern is the nitrogen.
        n_idx = match[0]
        n_atom = mol.GetAtomWithIdx(n_idx)
        
        # Check that the nitrogen has exactly two heavy-atom neighbors.
        if n_atom.GetDegree() != 2:
            continue  # not the expected substitution
        
        # Check that the nitrogen carries exactly one hydrogen.
        if n_atom.GetTotalNumHs() != 1:
            continue
        
        # Identify the backbone neighbor: in our SMARTS, the second atom (index 1)
        backbone_neighbor_idx = match[1]
        
        acyl_candidate = None
        # Look among the nitrogen bonds for the non-backbone neighbor.
        for bond in n_atom.GetBonds():
            nbr = bond.GetOtherAtom(n_atom)
            if nbr.GetIdx() == backbone_neighbor_idx:
                continue  # skip the backbone partner
            # The bond from the nitrogen to the acyl candidate must be a single bond.
            if bond.GetBondType() != Chem.BondType.SINGLE:
                continue
            acyl_candidate = nbr
            break
            
        if acyl_candidate is None:
            continue  # no valid acyl candidate found
        
        # The acyl candidate must be a carbon atom.
        if acyl_candidate.GetAtomicNum() != 6:
            continue
        
        # Check that the acyl candidate carbon is sp2-hybridized (typical for an amide carbonyl).
        if acyl_candidate.GetHybridization() != Chem.HybridizationType.SP2:
            continue
        
        # Check that acyl candidate carbon is double-bonded to exactly one oxygen.
        carbonyl_count = 0
        for bond in acyl_candidate.GetBonds():
            # Skip the bond to the nitrogen.
            other = bond.GetOtherAtom(acyl_candidate)
            if other.GetIdx() == n_idx:
                continue
            if bond.GetBondType() == Chem.BondType.DOUBLE and other.GetAtomicNum() == 8:
                carbonyl_count += 1
        if carbonyl_count != 1:
            continue

        # If all conditions are met, we have a sphinganine backbone with proper N-acyl substitution.
        return True, "Molecule contains sphinganine backbone with proper N-acyl substitution"
    
    return False, ("Sphinganine backbone found but nitrogen is not properly acylated "
                   "(substitution pattern does not match criteria)")

# For quick testing when running as a script.
if __name__ == "__main__":
    # Test one representative true positive example.
    test_smiles = "CCCCCCCCCCCCC[C@@H](O)[C@H](CO)NC(=O)CCCCCCCCC\\C=C/CCCCCCCC"
    result, reason = is_N_acylsphinganine(test_smiles)
    print("Test molecule:", result, reason)