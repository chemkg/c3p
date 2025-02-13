"""
Classifies: CHEBI:2440 acrovestone
"""
"""
Classifies: Acrovestone – A polyphenol isolated from Acronychia pedunculata that exhibits moderate antioxidant and antityrosinase activities.
Heuristic criteria used:
  - Presence of an isoflavone-like core (a benzopyranone fused to an aromatic ring)
  - At least two phenolic (aromatic –OH) groups
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_acrovestone(smiles: str):
    """
    Determines if a molecule is an acrovestone-type polyphenol based on its SMILES string.
    Heuristics:
      1. The molecule must contain an isoflavone-like core. A revised, more permissive SMARTS pattern is used
         to capture a benzopyranone system fused with an aromatic ring. This pattern should tolerate additional
         substituents such as sugars.
      2. The molecule should contain at least two phenolic hydroxyl groups (-OH attached to an aromatic carbon).

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule meets the heuristic criteria for acrovestone, False otherwise.
        str: A reason explaining the classification decision.
    """
    # Parse the input SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Remove explicit hydrogens to simplify substructure matching.
    mol = Chem.RemoveHs(mol)
    
    # Revised SMARTS for an isoflavone-like core:
    # This pattern is intended to capture a benzopyranone fused to an aromatic ring.
    # The pattern "c1ccc2oc(=O)cc2c1" allows for extra substituents in the test molecules to be tolerated.
    isoflavone_core_smarts = "c1ccc2oc(=O)cc2c1"
    core_query = Chem.MolFromSmarts(isoflavone_core_smarts)
    if core_query is None:
        return False, "Error in isoflavone core SMARTS pattern"
    
    if not mol.HasSubstructMatch(core_query):
        return False, "Isoflavone-like core not detected"
    
    # Define a SMARTS pattern for a phenolic hydroxyl group.
    # This pattern looks for an -OH group bonded to any aromatic carbon.
    phenol_smarts = "[$([c][OX2H])]"
    phenol_query = Chem.MolFromSmarts(phenol_smarts)
    if phenol_query is None:
        return False, "Error in phenol SMARTS pattern"
    
    # Count the number of phenolic hydroxyl groups found.
    phenol_matches = mol.GetSubstructMatches(phenol_query)
    if len(phenol_matches) < 2:
        return False, f"Found {len(phenol_matches)} phenolic -OH groups; at least 2 required"
    
    return True, "Contains an isoflavone-like core and sufficient phenolic hydroxyl groups indicative of acrovestone polyphenols"

# Example usage (for testing purposes)
if __name__ == "__main__":
    # Test one example from the provided list:
    test_smiles = "O1C(C(O)C(O)C(O)C1OC2=C(CCC(O)(C)C)C=3OC=C(C(=O)C3C(O)=C2)C4=CC=C(O)C=C4"
    result, reason = is_acrovestone(test_smiles)
    print("Acrovestone classification:", result)
    print("Reason:", reason)