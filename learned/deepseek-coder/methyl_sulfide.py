"""
Classifies: CHEBI:86315 methyl sulfide
"""
"""
Classifies: CHEBI:48375 methyl sulfide
"""
from rdkit import Chem

def is_methyl_sulfide(smiles: str):
    """
    Determines if a molecule is a methyl sulfide based on its SMILES string.
    A methyl sulfide is any aliphatic sulfide where at least one of the organyl groups attached to the sulfur is a methyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a methyl sulfide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the methyl sulfide pattern: sulfur bonded to a methyl group
    methyl_sulfide_pattern = Chem.MolFromSmarts("[SX2][CH3]")  # S bonded to a methyl group

    # Check if the molecule contains the methyl sulfide pattern
    matches = mol.GetSubstructMatches(methyl_sulfide_pattern)
    if not matches:
        return False, "No sulfur atom bonded to a methyl group found"

    # Ensure the sulfur is aliphatic (not in a ring or aromatic)
    for match in matches:
        sulfur_idx = match[0]  # Index of the sulfur atom in the match
        sulfur_atom = mol.GetAtomWithIdx(sulfur_idx)
        if not sulfur_atom.IsInRing() and not sulfur_atom.GetIsAromatic():
            return True, "Contains an aliphatic sulfur atom bonded to at least one methyl group"

    return False, "Sulfur atom is not aliphatic (may be in a ring or aromatic)"