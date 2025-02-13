"""
Classifies: CHEBI:46640 diketone
"""
"""
Classifies: A compound that contains exactly two ketone functionalities (diketone)
"""

from rdkit import Chem

def is_diketone(smiles: str):
    """
    Determines if a molecule is a diketone based on its SMILES string.
    A diketone is defined as a compound containing exactly two ketone groups.
    A ketone functionality is identified as a carbonyl (C=O) where the carbon
    is bound to two carbon atoms (i.e. excluding aldehydes).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a diketone, False otherwise
        str: Reason for classification or rejection
    """
    # Parse SMILES string into a molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a SMARTS pattern for a ketone group.
    # This pattern matches a carbonyl where the carbon (C ) is bound to two carbon atoms.
    ketone_pattern = Chem.MolFromSmarts("[#6][C](=O)[#6]")
    if ketone_pattern is None:
        return False, "Error in ketone SMARTS pattern"

    # Find all unique substructure matches of the ketone pattern in the molecule
    ketone_matches = mol.GetSubstructMatches(ketone_pattern)
    n_ketones = len(ketone_matches)

    # Check if there are exactly two ketone functionalities
    if n_ketones == 2:
        return True, "Compound contains exactly 2 ketone groups."
    else:
        return False, f"Compound contains {n_ketones} ketone group(s), which does not equal 2."

# Example usage (you can remove these lines if integrating into another codebase):
if __name__ == "__main__":
    test_smiles = [
        "O=C(CCCCCCCCCCCCCCC)CC(=O)CCCCCCCCC",  # 10,12-Heptacosanedione
        "CCC(=O)C(C)=O",                       # pentane-2,3-dione
        "O=C(CC(=O)c1ccccc1)c1ccccc1"           # dibenzoylmethane
    ]
    for smi in test_smiles:
        result, reason = is_diketone(smi)
        print(f"SMILES: {smi}")
        print(f"Result: {result}, Reason: {reason}")
        print("-"*50)