"""
Classifies: CHEBI:88061 polyamine
"""
"""
Classifies: Polyamine
Definition: Any organic amino compound that contains two or more amino groups.

Amino groups are identified by counting nitrogen atoms that have at least one hydrogen attached.
"""
from rdkit import Chem

def is_polyamine(smiles: str):
    """
    Determines if a molecule is a polyamine based on its SMILES string.
    A polyamine is defined as any organic amino compound that contains two or more amino groups.
    An "amino group" here is identified as a nitrogen (atomic number 7) that has at least one hydrogen attached.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polyamine, False otherwise.
        str: Reason for classification.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check that the molecule is organic (contains at least one carbon atom).
    if not any(atom.GetAtomicNum() == 6 for atom in mol.GetAtoms()):
        return False, "Molecule does not contain carbon atoms, thus not organic"

    # Add hydrogens so we have explicit counts of attached hydrogens.
    mol = Chem.AddHs(mol)

    # Count amino groups. We define an amino group as a nitrogen atom with at least one attached hydrogen.
    amino_count = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7:  # Only consider nitrogen atoms.
            # GetTotalNumHs returns the number of hydrogens attached (explicit + implicit).
            if atom.GetTotalNumHs() >= 1:
                amino_count += 1

    # Require at least two amino groups.
    if amino_count < 2:
        return False, f"Contains {amino_count} amino group(s), need at least 2 for polyamine"
    else:
        return True, f"Contains {amino_count} amino groups, satisfying polyamine criteria"