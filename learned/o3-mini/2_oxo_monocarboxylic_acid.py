"""
Classifies: CHEBI:35910 2-oxo monocarboxylic acid
"""
"""
Classifies: 2-oxo monocarboxylic acid
Definition: Any monocarboxylic acid having a 2-oxo substituent.
"""

from rdkit import Chem

def is_2_oxo_monocarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a 2-oxo monocarboxylic acid based on its SMILES string.
    This requires:
      1. The molecule to have exactly one carboxylic acid group (-C(=O)O[H]).
      2. The carbon attached to the carboxyl carbon (the alpha carbon) must contain a C=O bond,
         representing the "2-oxo" substituent.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a 2-oxo monocarboxylic acid, False otherwise.
        str: Reason for the classification decision.
    """
    # Parse the SMILES string into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a SMARTS pattern for a carboxylic acid group.
    # This pattern matches a carbon (sp2) doubly bonded to oxygen and single bonded to an -OH group.
    acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H]")
    acid_matches = mol.GetSubstructMatches(acid_pattern)
    
    # Ensure that exactly one carboxylic acid group is present (monocarboxylic acid)
    if len(acid_matches) == 0:
        return False, "No carboxylic acid group found"
    if len(acid_matches) > 1:
        return False, "More than one carboxylic acid group found; not a monocarboxylic acid"

    # The match returns a tuple of atom indices matching the pattern.
    # According to our SMARTS, index 0 is the carboxyl carbon.
    acid_carbon_idx = acid_matches[0][0]
    acid_carbon = mol.GetAtomWithIdx(acid_carbon_idx)

    # Identify the alpha carbon: it should be a neighbor of the carboxyl carbon that is not oxygen.
    alpha_atom = None
    for neighbor in acid_carbon.GetNeighbors():
        # Skip if the neighbor is oxygen (carboxyl oxygens)
        if neighbor.GetAtomicNum() == 8:
            continue
        alpha_atom = neighbor
        break
    if alpha_atom is None:
        return False, "No suitable alpha carbon found attached to the carboxyl group"

    # Check that the alpha carbon has a carbonyl substituent (C=O).
    # This is done by looking through the bonds of the alpha carbon to see if it is double bonded to an oxygen.
    has_oxo = False
    for bond in alpha_atom.GetBonds():
        # A double bond to oxygen qualifies as the 2-oxo group
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
            # Identify the neighboring atom in this bond
            other_atom = bond.GetOtherAtom(alpha_atom)
            if other_atom.GetAtomicNum() == 8:
                has_oxo = True
                break

    if not has_oxo:
        return False, "No 2-oxo substituent found on the alpha carbon of the carboxyl group"

    return True, "Found a monocarboxylic acid with a 2-oxo substituent on the alpha carbon"