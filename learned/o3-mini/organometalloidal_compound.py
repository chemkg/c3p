"""
Classifies: CHEBI:143084 organometalloidal compound
"""
"""
Classifies: Organometalloidal Compound
Definition: A compound having bonds between one or more metalloid atoms and one or more carbon atoms of an organyl group.
"""
from rdkit import Chem

def is_organometalloidal_compound(smiles: str):
    """
    Determines if a molecule is an organometalloidal compound based on its SMILES string.
    A compound is organometalloidal if it contains at least one bond between a metalloid atom and a carbon atom
    (typically present in an organyl group).
    
    For our implementation, metalloids are considered to be:
       Boron (B, atomic number 5),
       Silicon (Si, 14),
       Germanium (Ge, 32),
       Arsenic (As, 33),
       Antimony (Sb, 51),
       Tellurium (Te, 52),
       Polonium (Po, 84).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule qualifies as an organometalloidal compound, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the set of metalloids by their atomic numbers.
    metalloids = {5, 14, 32, 33, 51, 52, 84}
    
    # Iterate over all bonds in the molecule.
    for bond in mol.GetBonds():
        atom1 = bond.GetBeginAtom()
        atom2 = bond.GetEndAtom()
        atomic_num1 = atom1.GetAtomicNum()
        atomic_num2 = atom2.GetAtomicNum()
        
        # Check if one atom is a metalloid and the other is carbon (atomic number 6).
        if (atomic_num1 in metalloids and atomic_num2 == 6) or (atomic_num2 in metalloids and atomic_num1 == 6):
            return True, f"Found a bond between {atom1.GetSymbol()} and {atom2.GetSymbol()}"
    
    return False, "No bond found between a metalloid atom and a carbon atom of an organyl group."

# Example usage:
if __name__ == "__main__":
    # Test with methylarsonic acid: SMILES: "C[As](O)(O)=O"
    result, reason = is_organometalloidal_compound("C[As](O)(O)=O")
    print(result, reason)