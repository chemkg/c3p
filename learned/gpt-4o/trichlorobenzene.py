"""
Classifies: CHEBI:27096 trichlorobenzene
"""
from rdkit import Chem

def is_trichlorobenzene(smiles: str):
    """
    Determines if a molecule is a trichlorobenzene based on its SMILES string.
    A trichlorobenzene is any benzene ring with exactly three chlorine substituents.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule contains a trichlorobenzene, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define benzene ring pattern
    benzene_pattern = Chem.MolFromSmarts("c1ccccc1")

    # Find all benzene rings in the molecule
    benzene_matches = mol.GetSubstructMatches(benzene_pattern)

    if not benzene_matches:
        return False, "No benzene ring found"

    # Go through each benzene ring individually
    for match in benzene_matches:
        chlorine_count = 0  # Count of chlorine atoms connected directly to this benzene ring

        for atom_idx in match:
            atom = mol.GetAtomWithIdx(atom_idx)
            for neighbor in atom.GetNeighbors():
                if neighbor.GetSymbol() == 'Cl':
                    # Ensure chlorine is directly attached to this benzene and not involved in other rings
                    chlorine_count += 1

        # Check if exactly three chlorines are attached to this benzene ring
        if chlorine_count == 3:
            return True, "Contains a benzene ring with exactly three chlorine substituents"

    return False, "No benzene ring with exactly three chlorine substituents found"