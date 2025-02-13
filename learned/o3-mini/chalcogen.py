"""
Classifies: CHEBI:33303 chalcogen
"""
"""
Classifies: Any p-block element belonging to group 16 (chalcogens)
Group 16 elements: oxygen (8), sulfur (16), selenium (34), tellurium (52), polonium (84)
A valid chalcogen entity for this classifier must be a single free atom, without any bonds.
"""

from rdkit import Chem

def is_chalcogen(smiles: str):
    """
    Determines if a chemical entity, based on its SMILES string, is a free chalcogen element.
    A valid free chalcogen should be represented as a single atom SMILES (e.g., [Te], [S], [18O], etc.).
    
    Args:
        smiles (str): SMILES string of the chemical entity.
        
    Returns:
        bool: True if the entity is a free chalcogen element, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Ensure that the SMILES represents exactly one atom with no bonds.
    # Free atoms should have exactly one atom and no bonds.
    if mol.GetNumAtoms() != 1 or mol.GetNumBonds() > 0:
        return False, "SMILES does not represent a free single atom; it appears to be a molecule or ion."

    # Get the single atom in the molecule.
    atom = mol.GetAtomWithIdx(0)
    atomic_num = atom.GetAtomicNum()

    # Define the set of group 16 (chalcogen) atomic numbers
    chalcogens = {8, 16, 34, 52, 84}
    if atomic_num in chalcogens:
        return True, f"Element with atomic number {atomic_num} is a chalcogen."
    else:
        return False, f"Element with atomic number {atomic_num} is not a chalcogen."

# Example usage (for testing, uncomment the following lines):
# print(is_chalcogen("[S]"))       # Expected: True, because sulfur (atomic number 16) is a chalcogen.
# print(is_chalcogen("[H][Te-][H]")) # Expected: False, because it represents a compound not a free atom.