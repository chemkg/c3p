"""
Classifies: CHEBI:53156 polychlorobiphenyl
"""
"""
Classifies: Polychlorobiphenyl
Definition: A biphenyl compound containing between 2 and 10 chlorine atoms attached to the two benzene rings.
"""
from rdkit import Chem

def is_polychlorobiphenyl(smiles: str):
    """
    Determines if a molecule is a polychlorobiphenyl based on its SMILES string.
    A polychlorobiphenyl is defined as a biphenyl compound where between 2 and 10 chlorine atoms
    are directly attached to the two benzene rings.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is classified as a polychlorobiphenyl, False otherwise.
        str: Explanation for the classification decision.
    """
    
    # Parse the SMILES string into a molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a SMARTS pattern for a biphenyl scaffold: two benzene rings connected by a single bond.
    biphenyl_pattern = Chem.MolFromSmarts("c1ccccc1-c2ccccc2")
    if not mol.HasSubstructMatch(biphenyl_pattern):
        return False, "Molecule does not contain a biphenyl scaffold (two benzene rings connected by a single bond)."

    # Get the first biphenyl match (12 atoms: 6 from each benzene ring).
    matches = mol.GetSubstructMatches(biphenyl_pattern)
    if not matches:
        return False, "No valid biphenyl scaffold match found."

    # Use the first match found.
    biphenyl_atom_indices = set(matches[0])
    
    # Count chlorine atoms (atomic number 17) directly attached to the biphenyl rings.
    # Use a set to avoid double counting a chlorine attached to two ring carbons.
    cl_atoms_found = set()
    for idx in biphenyl_atom_indices:
        atom = mol.GetAtomWithIdx(idx)
        for neighbor in atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 17:
                cl_atoms_found.add(neighbor.GetIdx())

    cl_count = len(cl_atoms_found)
    
    # Check that the number of chlorine atoms falls between 2 and 10.
    if cl_count < 2:
        return False, f"Found {cl_count} chlorine substituent(s) on the biphenyl scaffold; at least 2 are required."
    if cl_count > 10:
        return False, f"Found {cl_count} chlorine substituent(s) on the biphenyl scaffold; no more than 10 are allowed."
    
    return True, f"Molecule is a polychlorobiphenyl with {cl_count} chlorine substituent(s) on its biphenyl scaffold."

# Example usage:
if __name__ == "__main__":
    test_smiles = "Clc1cc(Cl)c(Cl)c(c1)-c1cc(Cl)c(Cl)c1Cl"  # For 2,2',3,3',5,5'-hexachlorobiphenyl
    result, reason = is_polychlorobiphenyl(test_smiles)
    print(result, reason)