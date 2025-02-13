"""
Classifies: CHEBI:38958 indole alkaloid
"""
"""
Classifies: Indole alkaloid

Definition: An alkaloid containing an indole skeleton.
This implementation uses a SMARTS pattern to detect the indole core 
(c1ccc2[nH]ccc2c1) and then requires that at least one match corresponds
to a classical indole fragment (with 9 or 10 atoms in the match). 
It also includes a simple amide bond check to help filter out peptide‐like 
molecules (which sometimes contain a Trp residue but are not indole alkaloids).
"""

from rdkit import Chem

def is_indole_alkaloid(smiles: str):
    """
    Determines if a given molecule is an indole alkaloid based on its SMILES string.
    For our purposes, the molecule qualifies if it:
      (1) Can be parsed as a molecule,
      (2) Contains at least one nitrogen atom (as a rough check for alkaloid nature),
      (3) Contains an indole skeleton as defined by the SMARTS pattern "c1ccc2[nH]ccc2c1"
          with at least one match whose atom count is 9 (or at most one extra atom to allow slight extension),
      (4) Does NOT contain multiple amide bonds (which would suggest a peptide rather than an alkaloid).
      
    Args:
        smiles (str): The SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as an indole alkaloid, False otherwise.
        str: Explanation string for the decision.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check that the molecule has at least one nitrogen atom.
    num_nitrogen = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    if num_nitrogen == 0:
        return False, "Molecule does not contain any nitrogen atoms; unlikely to be an alkaloid"
    
    # (Optional) Filter out obvious peptides: count amide bonds via a simple SMARTS. 
    # Most indole alkaloids have at most one amide bond. (Peptides will contain two or more.)
    amide_smarts = Chem.MolFromSmarts("C(=O)N")
    amide_matches = mol.GetSubstructMatches(amide_smarts)
    if len(amide_matches) >= 2:
        return False, "Molecule contains multiple amide bonds; likely a peptide rather than an indole alkaloid"
    
    # Define the indole SMARTS pattern.
    indole_smarts = "c1ccc2[nH]ccc2c1"  # canonical indole substructure
    indole_query = Chem.MolFromSmarts(indole_smarts)
    if indole_query is None:
        return False, "Error in indole SMARTS pattern"
    
    # Look for indole substructure matches.
    matches = mol.GetSubstructMatches(indole_query)
    if not matches:
        return False, "No indole skeleton found"
    
    # Check that at least one match is of the expected size (9 atoms, possibly 10 to allow slight extensions).
    valid_match = False
    for match in matches:
        # For a classical indole, the match should yield 9 atoms.
        # However, sometimes one extra atom may be included due to fused extensions.
        if len(match) == 9 or len(match) == 10:
            valid_match = True
            break
    if not valid_match:
        return False, "Indole-like substructure found but does not match the expected indole fragment size"
    
    # Passed all tests: assume the molecule contains an indole skeleton and is of alkaloid type.
    return True, "Molecule contains an indole skeleton and minimal amide bonds, consistent with an indole alkaloid"

# For testing (if this code is run as a script), you can include a main section.
if __name__ == "__main__":
    # Test with one example (Ochropposinine)
    test_smiles = "CC[C@]1(CN2CCC=3C4=CC(=C(C=C4NC3[C@@]2(C[C@@]1(CCO)[H])[H])OC)OC)[H]"
    result, reason = is_indole_alkaloid(test_smiles)
    print("Result:", result)
    print("Reason:", reason)