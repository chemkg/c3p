"""
Classifies: CHEBI:33303 chalcogen
"""
"""
Classifies: Any p-block element belonging to group 16 (chalcogens)
Group 16 elements: O (8), S (16), Se (34), Te (52), Po (84)
"""

from rdkit import Chem

def is_chalcogen(smiles: str):
    """
    Determines if a molecule is a chalcogen element based on its SMILES string.
    A chalcogen is any p-block element belonging to group 16 of the periodic table:
    oxygen (O), sulfur (S), selenium (Se), tellurium (Te), or polonium (Po).

    Args:
        smiles (str): SMILES string of the chemical entity (typically an atom)

    Returns:
        bool: True if the entity is a chalcogen element, False otherwise.
        str: Description/reason for the classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Expect the molecule to be a single atom.
    if mol.GetNumAtoms() != 1:
        return False, "SMILES does not represent a single atom; it may be a molecule."

    # Get the single atom and its atomic number.
    atom = mol.GetAtomWithIdx(0)
    atomic_num = atom.GetAtomicNum()

    # Define atomic numbers for group 16 elements (chalcogens)
    chalcogens = {8, 16, 34, 52, 84}
    if atomic_num in chalcogens:
        return True, f"Element with atomic number {atomic_num} is a chalcogen."
    else:
        return False, f"Element with atomic number {atomic_num} is not a chalcogen."

# Example usage (uncomment for testing):
# print(is_chalcogen("[S]"))         # Expected: (True, "Element with atomic number 16 is a chalcogen.")
# print(is_chalcogen("[C]"))         # Expected: (False, "Element with atomic number 6 is not a chalcogen.")