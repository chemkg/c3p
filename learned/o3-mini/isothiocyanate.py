"""
Classifies: CHEBI:52221 isothiocyanate
"""
"""
Classifies: Isothiocyanate (organosulfur compound with the general formula R-N=C=S)
"""

from rdkit import Chem

def is_isothiocyanate(smiles: str):
    """
    Determines if a molecule is an isothiocyanate based on its SMILES string.
    An isothiocyanate is defined as an organosulfur compound with the general formula R-N=C=S.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is an isothiocyanate, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a SMARTS pattern for the isothiocyanate group.
    # This pattern requires a nitrogen double-bonded to a carbon,
    # which is double-bonded to a sulfur: N=C=S.
    iso_smarts = "N=C=S"
    iso_group = Chem.MolFromSmarts(iso_smarts)
    
    # If the molecule does not contain this substructure, it is not an isothiocyanate.
    if not mol.HasSubstructMatch(iso_group):
        return False, "No isothiocyanate functional group (N=C=S) found"

    # Get all matches for the isothiocyanate substructure.
    matches = mol.GetSubstructMatches(iso_group)
  
    # Check for at least one occurrence where the nitrogen atom has an additional substituent.
    # This is to confirm the "R" group in R-N=C=S. The nitrogen should have other connection(s)
    # besides the carbon participating in N=C=S.
    valid = False
    for match in matches:
        # match: tuple of indices (n_index, c_index, s_index) corresponding to N, C, S.
        n_atom = mol.GetAtomWithIdx(match[0])
        # Exclude the carbon that is part of the N=C=S group.
        neighbors = [nbr for nbr in n_atom.GetNeighbors() if nbr.GetIdx() != match[1]]
        if neighbors:  # Has at least one additional substituent
            valid = True
            break

    if not valid:
        return False, "Isothiocyanate group found but nitrogen lacks an R substituent"

    return True, "Contains isothiocyanate functional group (R-N=C=S) with appropriate substitution"