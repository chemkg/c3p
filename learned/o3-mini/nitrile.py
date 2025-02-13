"""
Classifies: CHEBI:18379 nitrile
"""
"""
Classifies: nitrile compounds (RC#N)
A nitrile is defined as a compound featuring a carbon–nitrogen triple bond 
with the carbon being substituted (i.e. not just HC#N).
"""

from rdkit import Chem

def is_nitrile(smiles: str):
    """
    Determines if a molecule is a nitrile based on its SMILES string.
    A nitrile should contain at least one C#N functional group where the C is substituted.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if molecule is a nitrile (contains at least one substituted nitrile group), False otherwise.
        str: Reason for classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the nitrile functional group pattern (a carbon triple-bonded to nitrogen)
    nitrile_pattern = Chem.MolFromSmarts("[C]#[N]")
    nitrile_matches = mol.GetSubstructMatches(nitrile_pattern)
    if not nitrile_matches:
        return False, "No nitrile (C#N) functional group found"

    # For each nitrile match, ensure the carbon atom is substituted (i.e., has a neighbor other than the nitrile nitrogen)
    for match in nitrile_matches:
        c_idx = match[0]  # index of carbon in the nitrile group
        nitrile_c = mol.GetAtomWithIdx(c_idx)

        # Count explicit neighbors; note that RDKit does not include implicit hydrogens.
        # In a pure HCN fragment, the only explicit neighbor of the nitrile carbon is the nitrile nitrogen.
        # For a substituted nitrile (RC#N), the nitrile carbon should have at least one other neighbor.
        neighbor_atoms = nitrile_c.GetNeighbors()
        # Count neighbors that are not the nitrile nitrogen (atomic num 7) because that is part of the nitrile group.
        substituent_neighbors = [atom for atom in neighbor_atoms if atom.GetAtomicNum() != 7]

        if substituent_neighbors:
            # Found a nitrile group where the carbon is substituted by a non-hydrogen atom.
            return True, "Molecule contains a substituted nitrile group (RC#N)"

    # If we reach this point, then either only HCN groups were found or the nitrile carbon is not substituted.
    return False, "Only unsubstituted nitrile group(s) (HC#N) found; no C-substituted nitrile group present"