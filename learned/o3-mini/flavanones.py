"""
Classifies: CHEBI:28863 flavanones
"""
"""
Classifies: Flavanones
Definition: Members of the class of flavans with a 3,4-dihydro-2-aryl-2H-1-benzopyran-4-one skeleton and its substituted derivatives.
Improvement:
  - Uses an improved SMARTS pattern that forces the two ring carbons immediately following the carbonyl to be in a ring and non‐aromatic.
  - Checks that these two carbons are sp3 hybridized (thus confirming the “dihydro” character).
  - Uses useChirality=False in substructure matching.
"""

from rdkit import Chem
from rdkit.Chem import rdchem

def is_flavanones(smiles: str):
    """
    Determines if a molecule is a flavanone (or a substituted derivative) based on its SMILES string.
    A flavanone must contain the 3,4-dihydro-2-aryl-2H-1-benzopyran-4-one core.

    The approach is as follows:
      1. Parse the SMILES into an RDKit Mol object.
      2. Define a SMARTS pattern for the flavanone core.
         The pattern is "O=C1[C;R;!a][C;R;!a](c2ccccc2)Oc2ccccc12" which means:
           • "O=C1"  : a carbonyl group (the C4=O)
           • "[C;R;!a]" : a ring carbon that is not aromatic (the saturated carbon at position 3)
           • "[C;R;!a](c2ccccc2)" : a second ring carbon (position 2) that carries an aryl substituent,
           • "Oc2ccccc12" : an oxygen bridging to another aryl ring (the fused benzene ring).
      3. The matching is done ignoring chirality by setting useChirality=False.
      4. For every match the two saturated carbons (positions 1 and 2 in the SMARTS, after the carbonyl) are checked.
         They must be sp3 hybridized to confirm the “dihydro” nature.
         
    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule contains the flavanone core (3,4-dihydro-2-aryl-2H-1-benzopyran-4-one); 
              False otherwise.
        str: Explanation of the classification.
    """
    
    # Parse SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define an improved SMARTS pattern:
    # Instead of [C;!R] (which means not in a ring), we require the atom to be in a ring and non-aromatic: [C;R;!a]
    flavanone_pattern = Chem.MolFromSmarts("O=C1[C;R;!a][C;R;!a](c2ccccc2)Oc2ccccc12")
    if flavanone_pattern is None:
        return False, "Failed to create SMARTS pattern for flavanone core"
    
    # Find all matches in the molecule while ignoring chirality
    matches = mol.GetSubstructMatches(flavanone_pattern, useChirality=False)
    if not matches:
        return False, "Core flavanone skeleton (3,4-dihydro-2-aryl-2H-1-benzopyran-4-one) not found"
    
    # For each match, check that the two saturated carbons are indeed sp3 hybridized.
    # Our SMARTS is structured as follows:
    # match[0]: carbonyl carbon (expected to be sp2)
    # match[1]: the first saturated carbon (should be sp3)
    # match[2]: the second saturated carbon (should be sp3, bearing the aryl substituent)
    valid_match_found = False
    for match in matches:
        if len(match) < 3:
            continue  # safety check
        
        atom1 = mol.GetAtomWithIdx(match[1])
        atom2 = mol.GetAtomWithIdx(match[2])
        if (atom1.GetHybridization() == rdchem.HybridizationType.SP3 and
            atom2.GetHybridization() == rdchem.HybridizationType.SP3):
            valid_match_found = True
            break
    if not valid_match_found:
        return False, "Substructure match found but one or both saturated atoms are not sp3; likely not a true flavanone core"
    
    return True, "Molecule contains the flavanone core skeleton (3,4-dihydro-2-aryl-2H-1-benzopyran-4-one)"

# Example usage and testing:
if __name__ == "__main__":
    # Test with a known flavanone example: naringenin
    test_smiles = "Oc1ccc(cc1)C1CC(=O)c2c(O)cc(O)cc2O1"
    result, explanation = is_flavanones(test_smiles)
    print(f"SMILES: {test_smiles}\nClassification: {result}\nReason: {explanation}")