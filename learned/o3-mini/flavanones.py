"""
Classifies: CHEBI:28863 flavanones
"""
"""
Classifies: Flavanones
Definition: Members of the class of flavans with a 3,4-dihydro-2-aryl-2H-1-benzopyran-4-one skeleton and its substituted derivatives.
Improvement:
  - Uses the known SMARTS pattern for the 2-phenylchroman-4-one (flavanone) core.
  - Specifies useChirality=False in substructure matching.
  - Inspects the match to ensure that the two saturated carbons in the heterocycle remain non-aromatic,
    as required for the dihydro (flavanone) core.
"""

from rdkit import Chem

def is_flavanones(smiles: str):
    """
    Determines if a molecule is a flavanone (or a substituted derivative) based on its SMILES string.
    A flavanone must contain the 3,4-dihydro-2-aryl-2H-1-benzopyran-4-one core.
    
    To improve upon a simple substructure match, this function:
      1. Parses the molecule and creates a SMARTS pattern based on the known flavanone core.
         (Pattern: "O=C1CC(c2ccccc2)Oc2ccccc12")
      2. Uses useChirality=False to allow variant stereochemistry.
      3. For every match, checks that the two carbons in the heterocycle (immediately following the carbonyl)
         are not aromatic (as they are sp3 in a true flavanone).
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a flavanone; False otherwise.
        str: Explanation of the reasoning.
    """
    
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern for the flavanone core:
    # This pattern attempts to capture the 3,4-dihydro-2-aryl-2H-1-benzopyran-4-one scaffold.
    # Note: The pattern "O=C1CC(c2ccccc2)Oc2ccccc12" was chosen because it represents a carbonyl (O=C1),
    # two connected saturated carbons (CC) with an aryl substitution on the second one,
    # and an oxygen bridging into a fused aromatic ring.
    flavanone_pattern = Chem.MolFromSmarts("O=C1CC(c2ccccc2)Oc2ccccc12")
    if flavanone_pattern is None:
        return False, "Failed to create SMARTS pattern for flavanone core"
    
    # Look for substructure matches ignoring chirality (to allow substituted derivatives)
    matches = mol.GetSubstructMatches(flavanone_pattern, useChirality=False)
    if not matches:
        return False, "Core flavanone skeleton (3,4-dihydro-2-aryl-2H-1-benzopyran-4-one) not found"
    
    # For each match, check that the two saturated carbons (the 1st and 2nd carbons in the heterocycle)
    # are indeed non-aromatic.
    # The SMARTS "O=C1CC(c2ccccc2)Oc2ccccc12" is interpreted by RDKit in the following approximate order:
    #   match[0]: the carbonyl carbon (sp2, may be aromatic–it is part of a carbonyl)
    #   match[1]: first ring carbon (should be saturated, non-aromatic)
    #   match[2]: second ring carbon (should be saturated, non-aromatic; also bears the aryl substituent)
    # The remaining atoms belong to aromatic rings.
    valid_match_found = False
    for match in matches:
        if len(match) < 3:
            # Sanity; our pattern should give at least 3 atoms in order.
            continue
        atom1 = mol.GetAtomWithIdx(match[1])
        atom2 = mol.GetAtomWithIdx(match[2])
        # In a true flavanone, these two carbons (positions 1 and 2 in our match) must be non-aromatic.
        if not atom1.GetIsAromatic() and not atom2.GetIsAromatic():
            valid_match_found = True
            break
    
    if not valid_match_found:
        return False, "Substructure match found but core saturated atoms are aromatic; likely not a true flavanone core"
    
    return True, "Molecule contains the flavanone core skeleton (3,4-dihydro-2-aryl-2H-1-benzopyran-4-one)"

# Example usage:
if __name__ == "__main__":
    # one of the provided examples: naringenin
    test_smiles = "Oc1ccc(cc1)C1CC(=O)c2c(O)cc(O)cc2O1"
    is_flava, reason = is_flavanones(test_smiles)
    print(f"SMILES: {test_smiles}\nClassification: {is_flava}\nReason: {reason}")