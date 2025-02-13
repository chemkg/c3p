"""
Classifies: CHEBI:2440 acrovestone
"""
"""
Classifies: Acrovestone – A polyphenol isolated from Acronychia pedunculata that exhibits moderate antioxidant and antityrosinase activities.
Heuristic criteria used:
  1. A more permissive detection of a benzopyranone-like (isoflavone/flavone-like) core.
  2. The presence of at least two free phenolic (aromatic –OH) groups.
Note: This heuristic is intended to work with decorated structures (e.g. glycosylated) that retain the characteristic core.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_acrovestone(smiles: str):
    """
    Determines if a molecule belongs to the acrovestone-like polyphenols based on its SMILES string.
    Heuristics:
      1. The molecule must contain a benzopyranone-like core (a flavone-like core) as captured by a more permissive SMARTS pattern.
      2. The molecule should contain at least two free phenolic hydroxyl groups (-OH attached directly to an aromatic carbon).

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

    # Optionally remove explicit hydrogens to simplify pattern matching.
    mol = Chem.RemoveHs(mol)

    # Revised SMARTS for a benzopyranone-like core (a deregulated flavone/isoflavone core).
    # This pattern is more permissive and should be able to match a benzopyranone even when substituents (e.g. sugars) are attached.
    core_smarts = "O=c1cc2oc(=O)cc2c1"
    core_query = Chem.MolFromSmarts(core_smarts)
    if core_query is None:
        return False, "Error in core SMARTS pattern"

    if not mol.HasSubstructMatch(core_query):
        return False, "Isoflavone-like core not detected"

    # Define a SMARTS pattern for a phenolic hydroxyl group.
    # This pattern finds an -OH group directly bound to an aromatic carbon.
    phenol_smarts = "[$([c][OX2H])]"   # matches an aromatic carbon attached to an OH.
    phenol_query = Chem.MolFromSmarts(phenol_smarts)
    if phenol_query is None:
        return False, "Error in phenol SMARTS pattern"

    phenol_matches = mol.GetSubstructMatches(phenol_query)
    if len(phenol_matches) < 2:
        return False, f"Found {len(phenol_matches)} phenolic -OH groups; at least 2 required"

    return True, "Contains a benzopyranone (isoflavone-like) core and sufficient free phenolic hydroxyl groups indicative of acrovestone polyphenols"

# Example usage (for testing purposes)
if __name__ == "__main__":
    # Test one example from the provided list:
    test_smiles = "O1C(C(O)C(O)C(O)C1OC2=C(CCC(O)(C)C)C=3OC=C(C(=O)C3C(O)=C2)C4=CC=C(O)C=C4"
    result, reason = is_acrovestone(test_smiles)
    print("Acrovestone classification:", result)
    print("Reason:", reason)