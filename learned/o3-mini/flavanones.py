"""
Classifies: CHEBI:28863 flavanones
"""
"""
Classifies: Flavanones
Definition: Members of the class of flavans with a 3,4-dihydro-2-aryl-2H-1-benzopyran-4-one skeleton and its substituted derivatives.
Improvement:
  - Uses an improved SMARTS pattern that forces the two ring carbons adjacent to the carbonyl to be non-aromatic.
  - Checks that these carbons have sp3 hybridization.
  - Uses useChirality=False in substructure matching.
"""

from rdkit import Chem
from rdkit.Chem import rdchem

def is_flavanones(smiles: str):
    """
    Determines if a molecule is a flavanone (or a substituted derivative) based on its SMILES string.
    A flavanone must contain the 3,4-dihydro-2-aryl-2H-1-benzopyran-4-one core.

    To improve upon a simple substructure match, this function:
      1. Parses the SMILES string and creates a SMARTS pattern based on the known flavanone core.
         The pattern used is "O=C1[C;!R][C;!R](c2ccccc2)Oc2ccccc12", where the [C;!R] tokens force the
         two carbons in the heterocycle immediately following the carbonyl to be non-aromatic.
      2. Uses useChirality=False to allow variant stereochemistry.
      3. After matching, inspects the qualified atoms to ensure they are sp3 hybridized.
         This further helps confirm that the two carbons are saturated (dihydro).
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a flavanone; False otherwise.
        str: Explanation of the reasoning.
    """
    
    # Attempt to parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define an improved SMARTS pattern for the flavanone core.
    # "O=C1[C;!R][C;!R](c2ccccc2)Oc2ccccc12"
    # Explanation:
    #   O=C1          -> A carbonyl attached to ring index 1.
    #   [C;!R]        -> a carbon that is forced to be non-ring aromatic (non-aromatic/saturated).
    #   [C;!R](c2ccccc2) -> a second saturated carbon with an aryl substituent.
    #   Oc2ccccc12   -> An oxygen bridging into a fused benzene ring.
    flavanone_pattern = Chem.MolFromSmarts("O=C1[C;!R][C;!R](c2ccccc2)Oc2ccccc12")
    if flavanone_pattern is None:
        return False, "Failed to create SMARTS pattern for flavanone core"
    
    # Find all substructure matches ignoring chirality.
    matches = mol.GetSubstructMatches(flavanone_pattern, useChirality=False)
    if not matches:
        return False, "Core flavanone skeleton (3,4-dihydro-2-aryl-2H-1-benzopyran-4-one) not found"
    
    # For each match, check that the two saturated carbons (which should be at positions 1 and 2 in our SMARTS)
    # are indeed sp3 hybridized. This ensures the "dihydro" character of the heterocycle.
    valid_match_found = False
    for match in matches:
        # Our SMARTS must return at least 3 atoms:
        # match[0] = carbonyl carbon,
        # match[1] = first (saturated) ring carbon,
        # match[2] = second (saturated) ring carbon (which bears the aryl substituent).
        if len(match) < 3:
            continue
        # Get the two candidate atoms.
        atom1 = mol.GetAtomWithIdx(match[1])
        atom2 = mol.GetAtomWithIdx(match[2])
        # Check hybridization: expected to be sp3
        if atom1.GetHybridization() != rdchem.HybridizationType.SP3 or atom2.GetHybridization() != rdchem.HybridizationType.SP3:
            continue  # Not a valid flavanone core if these carbons are not sp3.
        valid_match_found = True
        break

    if not valid_match_found:
        return False, "Substructure match found but one or both saturated atoms are not sp3; likely not a true flavanone core"
    
    return True, "Molecule contains the flavanone core skeleton (3,4-dihydro-2-aryl-2H-1-benzopyran-4-one)"

# Example usage:
if __name__ == "__main__":
    # Test with one of the provided examples: naringenin
    test_smiles = "Oc1ccc(cc1)C1CC(=O)c2c(O)cc(O)cc2O1"
    is_flava, reason = is_flavanones(test_smiles)
    print(f"SMILES: {test_smiles}\nClassification: {is_flava}\nReason: {reason}")