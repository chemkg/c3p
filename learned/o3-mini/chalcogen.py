"""
Classifies: CHEBI:33303 chalcogen
"""
"""
Classifies: Any p-block element belonging to group 16 (chalcogens)
Group 16 elements: oxygen (8), sulfur (16), selenium (34), tellurium (52), polonium (84)
A valid free chalcogen should be represented as a single atom with no bonds (after removing explicit hydrogens).
"""

from rdkit import Chem

def is_chalcogen(smiles: str):
    """
    Determines if a chemical entity is a free chalcogen element based on its SMILES string.
    We consider a valid free chalcogen representation to be one in which the heavy atom is a group 16 element,
    and it is not bonded to any other heavy atoms (after removing explicit hydrogens).

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
    
    # Remove explicit hydrogens so that we can check the connectivity of heavy atoms.
    mol_noHs = Chem.RemoveHs(mol)
    
    # A free atom (after removing hydrogens) must consist of exactly one atom and no bonds.
    if mol_noHs.GetNumAtoms() != 1 or mol_noHs.GetNumBonds() != 0:
        return False, "SMILES does not represent a free single atom (after removing hydrogens), appears to be a molecule or ion."
    
    # Get the single heavy atom.
    atom = mol_noHs.GetAtomWithIdx(0)
    atomic_num = atom.GetAtomicNum()
    
    # Define the set of group 16 (chalcogen) atomic numbers.
    chalcogens = {8, 16, 34, 52, 84}
    
    # Check if the heavy atom belongs to the chalcogens.
    if atomic_num in chalcogens:
        return True, f"Element with atomic number {atomic_num} is a chalcogen."
    else:
        return False, f"Element with atomic number {atomic_num} is not a chalcogen."

# Example usage (for testing, uncomment the following lines):
# print(is_chalcogen("[Te]"))         # Expected: True, tellurium is a chalcogen.
# print(is_chalcogen("[H][Te-][H]"))   # Expected: False, because after removing H atoms it would yield a free Te atom bound in a compound.