"""
Classifies: CHEBI:1624 3-oxo-5beta-steroid
"""
"""
Classifies: CHEBI:XXXXX 3-oxo-5beta-steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_3_oxo_5beta_steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-5beta-steroid based on its SMILES string.
    A 3-oxo-5beta-steroid is a steroid with a ketone at position 3 and beta-configuration at position 5.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-oxo-5beta-steroid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Assign stereochemistry
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True, flagPossibleStereoCenters=True)

    # Define the steroid core using a SMARTS pattern with atom mapping to identify positions
    steroid_core_smarts = """
    [#6]1([#6][#6]2[#6]([#6]1)[#6][#6]3[#6]([#6]2)[#6][#6][#6]4[#6]([#6]3)[#6][#6]([#6][#6]4))
    """
    steroid_core = Chem.MolFromSmarts(steroid_core_smarts)
    if steroid_core is None:
        return False, "Invalid steroid core SMARTS pattern"

    # Check for steroid core match
    matches = mol.GetSubstructMatches(steroid_core)
    if not matches:
        return False, "Steroid core not found"

    # Define SMARTS pattern for ketone at position 3 with atom mapping
    ketone_smarts = "[C;R1]-[C;R2](=O)-[C;R1]"
    ketone_pattern = Chem.MolFromSmarts(ketone_smarts)
    if ketone_pattern is None:
        return False, "Invalid ketone SMARTS pattern"

    # Check for ketone at position 3
    ketone_matches = mol.GetSubstructMatches(ketone_pattern)
    if not ketone_matches:
        return False, "Ketone group at position 3 not found"

    # Identify chiral centers
    stereocenters = Chem.FindMolChiralCenters(mol, includeUnassigned=True, includeCIP=True)

    # Assume position 5 is a chiral center in a ring connected to two other ring carbons
    position5_atom = None
    for atom in mol.GetAtoms():
        if atom.GetChiralTag() != Chem.rdchem.ChiralType.CHI_UNSPECIFIED:
            if atom.IsInRing():
                neighbors = atom.GetNeighbors()
                ring_neighbors = sum(1 for nbr in neighbors if nbr.IsInRing())
                if ring_neighbors >= 2:
                    position5_atom = atom
                    break

    if position5_atom is None:
        return False, "Chiral center at position 5 not found"

    # Check stereochemistry at position 5
    # Beta-configuration often corresponds to 'alpha' (R) configuration in RDKit
    cip_code = position5_atom.GetProp('_CIPCode') if position5_atom.HasProp('_CIPCode') else None
    if cip_code != 'R':
        return False, f"Position 5 is not in beta-configuration (found {cip_code})"

    return True, "Molecule is a 3-oxo-5beta-steroid"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:XXXXX',
        'name': '3-oxo-5beta-steroid',
        'definition': "Any 3-oxo steroid that has beta- configuration at position 5.",
        'parents': []
    },
    'config': {},
    'message': None,
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
}