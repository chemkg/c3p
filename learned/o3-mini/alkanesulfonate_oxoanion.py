"""
Classifies: CHEBI:134249 alkanesulfonate oxoanion
"""
"""
Classifies: alkanesulfonate oxoanion
Definition: 'An alkanesulfonate in which the carbon at position 1 is attached to R, which can represent hydrogens, a carbon chain, or other groups.'
This function looks for an sp3 carbon directly attached to a sulfonate moiety: S(=O)(=O)[O-].
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_alkanesulfonate_oxoanion(smiles: str):
    """
    Determines if a molecule contains an alkanesulfonate oxoanion moiety.
    This is defined as a sulfonate group –S(=O)(=O)[O-] attached to an sp3 carbon (i.e. an alkane fragment).

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if molecule contains the alkanesulfonate oxoanion feature, False otherwise.
        str: Reason for classification.
    """
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the SMARTS pattern:
    # [CX4] matches an sp3 carbon.
    # S(=O)(=O)[O-] matches a sulfonate group with two double-bonded oxygens and one negatively charged oxygen.
    alkanesulfonate_pattern = Chem.MolFromSmarts("[CX4]-S(=O)(=O)[O-]")
    if alkanesulfonate_pattern is None:
        return False, "Error constructing SMARTS pattern"

    # Search for the pattern in the molecule.
    matches = mol.GetSubstructMatches(alkanesulfonate_pattern)
    if not matches:
        return False, "No alkanesulfonate oxoanion substructure [CX4]-S(=O)(=O)[O-] found"
    
    # Optionally, further checks could be performed (for example, ensuring the carbon attached to S is not part of an aromatic ring,
    # or validating connectivity in complex molecules). For now, the presence of the SMARTS match is taken as sufficient.
    
    return True, f"Found {len(matches)} alkanesulfonate oxoanion instance(s) in the molecule"

# Below is a simple test harness that could be uncommented for testing purposes.
# if __name__ == "__main__":
#    test_smiles = [
#         "CS([O-])(=O)=O",  # methanesulfonate
#         "OCCN(CCO)CCS([O-])(=O)=O",  # 2-[bis(2-hydroxyethyl)amino]ethanesulfonate
#         "FC(S([O-])(=O)=O)(F)F",  # triflate
#         "CCCC",  # non-sulfonate example
#     ]
#     for smi in test_smiles:
#         result, reason = is_alkanesulfonate_oxoanion(smi)
#         print(f"SMILES: {smi} -> {result}, {reason}")