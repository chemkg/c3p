"""
Classifies: CHEBI:88061 polyamine
"""
"""
Classifies: Polyamine
Definition: Any organic amino compound that contains two or more amino groups.
An amino group is identified here as a nitrogen atom having at least one hydrogen atom directly attached.
"""
from rdkit import Chem

def is_polyamine(smiles: str):
    """
    Determines if a molecule is a polyamine based on its SMILES string.
    A polyamine is defined as any organic amino compound that contains two or more amino groups.
    An amino group is identified as a nitrogen (atomic number 7) that has at least one hydrogen atom attached.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polyamine, False otherwise.
        str: Reason for classification.
    """
    # Parse the SMILES string into a molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Ensure the molecule is organic (contains at least one carbon atom).
    if not any(atom.GetAtomicNum() == 6 for atom in mol.GetAtoms()):
        return False, "Molecule does not contain carbon atoms, thus not organic"
    
    # Make hydrogen atoms explicit.
    mol = Chem.AddHs(mol)
    
    # Count the number of amino groups by checking for nitrogen atoms that have at least one hydrogen neighbor.
    amino_count = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 7:
            continue
        # Count directly attached hydrogen atoms.
        h_neighbors = sum(1 for neighbor in atom.GetNeighbors() if neighbor.GetAtomicNum() == 1)
        if h_neighbors >= 1:
            amino_count += 1

    # At least two amino groups are required.
    if amino_count < 2:
        return False, f"Contains {amino_count} amino group(s), need at least 2 for polyamine"
    else:
        return True, f"Contains {amino_count} amino groups, satisfying polyamine criteria"