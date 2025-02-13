"""
Classifies: CHEBI:33303 chalcogen
"""
"""
Classifies: Any p-block element belonging to the group 16 family of the periodic table.
Group 16 elements: oxygen (8), sulfur (16), selenium (34), tellurium (52), polonium (84)
A valid free chalcogen should be represented as a single heavy atom with no bonds.
This improved version checks the original connectivity: after removing hydrogens,
the molecule must consist of exactly one atom and no bonds.
"""

from rdkit import Chem

def is_chalcogen(smiles: str):
    """
    Determines if a chemical entity is a free chalcogen element based on its SMILES string.
    A valid free chalcogen is defined as a single heavy atom (i.e. not hydrogen) with no bonds,
    ensuring that the representation is not part of a larger compound.

    Args:
        smiles (str): SMILES string of the chemical entity.
        
    Returns:
        bool: True if the entity is a free chalcogen element, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Remove explicit hydrogens to focus on the heavy atom connectivity.
    mol_noHs = Chem.RemoveHs(mol)
    
    # A free element, after hydrogen removal, must have exactly one atom and no bonds.
    if mol_noHs.GetNumAtoms() != 1:
        return False, "More than one heavy atom found after removing hydrogens; not a free element."
    
    if mol_noHs.GetNumBonds() != 0:
        return False, "Heavy atom shows bonding; not a free atom."
    
    # Get the single heavy atom.
    atom = mol_noHs.GetAtomWithIdx(0)
    atomic_num = atom.GetAtomicNum()
    
    # Check that the heavy atom is not hydrogen.
    if atomic_num == 1:
        return False, "The only atom present is hydrogen, not a chalcogen."
    
    # Define the set of chalcogen atomic numbers.
    chalcogens = {8, 16, 34, 52, 84}
    if atomic_num in chalcogens:
        return True, f"Element with atomic number {atomic_num} is a chalcogen."
    else:
        return False, f"Element with atomic number {atomic_num} is not a chalcogen."

# Example usage (for testing, uncomment the following lines):
# print(is_chalcogen("[Te]"))         # Expected: True, tellurium is a chalcogen.
# print(is_chalcogen("[H][Te-][H]"))   # Expected: False, because the SMILES represents a bonded species.