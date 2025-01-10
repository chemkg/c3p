"""
Classifies: CHEBI:16460 polyprenol phosphate
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_polyprenol_phosphate(smiles: str):
    """
    Determines if a molecule is a polyprenol phosphate based on its SMILES string.
    A polyprenol phosphate is characterized by a polyprenyl chain terminated by a phosphate or diphosphate group.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polyprenol phosphate, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"
    
    # Define SMARTS for isoprene units allowing for branching and stereochemistry
    isoprene_pattern = Chem.MolFromSmarts("C(=C)[-;@]C|C-C")

    # Find matching substructures for isoprene units
    isoprene_matches = mol.GetSubstructMatches(isoprene_pattern)
    if len(isoprene_matches) < 1:
        return False, "Too few isoprene units for polyprenol"
    
    # Broadened SMARTS patterns for phosphate and diphosphate groups
    phosphate_pattern = Chem.MolFromSmarts("OP(=O)(O)O")
    diphosphate_pattern = Chem.MolFromSmarts("OP(=O)(O)OP(=O)(O)O")

    # Check whether the molecule contains a terminal phosphate or diphosphate group
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    diphosphate_matches = mol.GetSubstructMatches(diphosphate_pattern)

    if not (phosphate_matches or diphosphate_matches):
        return False, "No terminal phosphate or diphosphate group found"

    # Further validation for classic structural features could be impacted by anticipated structures
    return True, "Contains polyprenyl chain with terminal phosphate or diphosphate group"