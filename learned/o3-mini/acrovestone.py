"""
Classifies: CHEBI:2440 acrovestone
"""
"""
Classifies: Acrovestone – A polyphenol isolated from Acronychia pedunculata that exhibits moderate antioxidant and antityrosinase activities.
Heuristic criteria used:
  - Presence of an isoflavone-like core (fused benzopyranone with an additional aromatic ring)
  - At least two phenolic (aromatic –OH) groups
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_acrovestone(smiles: str):
    """
    Determines if a molecule is an acrovestone-type polyphenol based on its SMILES string.
    Heuristics:
      1. The molecule must contain an isoflavone-like core. The SMARTS used here (c1ccc2c(c1)oc(=O)c3ccccc23)
         is meant to capture the benzopyranone scaffold fused to an aromatic ring.
      2. The molecule should have at least two phenolic hydroxyl groups (an OH group directly attached to an aromatic carbon).
         This often contributes to antioxidant activity.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule meets heuristic criteria for acrovestone, False otherwise.
        str: Reason for the classification decision.
    """
    # Try to parse the molecule from SMILES.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS for an isoflavone-like core.
    # This pattern represents a fused bicyclic aromatic ketone (benzopyran-4-one) connected to a second aromatic ring.
    # Note: This is a heuristic pattern.
    isoflavone_core_smarts = "c1ccc2c(c1)oc(=O)c3ccccc23"
    core_query = Chem.MolFromSmarts(isoflavone_core_smarts)
    if core_query is None:
        return False, "Error in SMARTS pattern"
    
    if not mol.HasSubstructMatch(core_query):
        return False, "Isoflavone-like core not detected"
    
    # Define a SMARTS pattern for a phenolic hydroxyl group.
    # This pattern looks for an -OH group attached to an aromatic carbon.
    phenol_smarts = "[cH][OX2H]"
    phenol_query = Chem.MolFromSmarts(phenol_smarts)
    if phenol_query is None:
        return False, "Error in phenol SMARTS pattern"
    
    phenol_matches = mol.GetSubstructMatches(phenol_query)
    if len(phenol_matches) < 2:
        return False, f"Found {len(phenol_matches)} phenolic -OH groups; at least 2 required"
    
    return True, "Contains isoflavone core with sufficient phenolic groups typical of acrovestone polyphenols"

# Example usage (for testing purposes)
if __name__ == "__main__":
    # (One among many examples provided)
    test_smiles = "O1C(C(O)C(O)C(O)C1OC2=C(CCC(O)(C)C)C=3OC=C(C(=O)C3C(O)=C2)C4=CC=C(O)C=C4"
    result, reason = is_acrovestone(test_smiles)
    print("Acrovestone classification:", result)
    print("Reason:", reason)