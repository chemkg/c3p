"""
Classifies: CHEBI:133004 bisbenzylisoquinoline alkaloid
"""
"""
Classifies: Bisbenzylisoquinoline alkaloid
Definition: A type of benzylisoquinoline alkaloid whose structures are built up of two benzylisoquinoline units 
linked by ether bridges. Various structural patterns resulting from additional bridging between the two units 
(by direct C–C linking or by methylenedioxy groups) are common.
"""
from rdkit import Chem

def is_bisbenzylisoquinoline_alkaloid(smiles: str):
    """
    Determines if a molecule is a bisbenzylisoquinoline alkaloid based on its SMILES string.
    The classification is made by checking two features:
      1. The molecule must contain at least two distinct benzylisoquinoline (tetrahydroisoquinoline)
         units. Here, we use the SMARTS "c1ccc2CCNCC2c1" as an approximate pattern.
      2. The two units must be connected by a bridging group. We check for an aromatic ether bridge,
         i.e. an oxygen atom with a bond to two aromatic carbons (SMARTS: "[c][O][c]").
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is classified as a bisbenzylisoquinoline alkaloid, False otherwise.
        str: Reason for classification.
    """
    # Parse the SMILES into a molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for a tetrahydroisoquinoline unit (approximation of a benzylisoquinoline unit)
    # This pattern represents an aromatic ring (benzene) fused to a saturated ring containing nitrogen.
    isoquinoline_smarts = "c1ccc2CCNCC2c1"
    isoquinoline_pattern = Chem.MolFromSmarts(isoquinoline_smarts)
    iso_matches = mol.GetSubstructMatches(isoquinoline_pattern)
    
    # Filter matches to obtain distinct (non‐overlapping) units.
    distinct_units = []
    for match in iso_matches:
        match_set = set(match)
        # Only add if this match does not share atoms with any unit already found.
        if not any(match_set & existing for existing in distinct_units):
            distinct_units.append(match_set)
    
    if len(distinct_units) < 2:
        return False, f"Found {len(distinct_units)} benzylisoquinoline unit(s); need at least 2"

    # Check for bridging: look for an aromatic ether bridge where an oxygen connects two aromatic carbons.
    ether_bridge_smarts = "[c][O][c]"
    ether_bridge = Chem.MolFromSmarts(ether_bridge_smarts)
    if not mol.HasSubstructMatch(ether_bridge):
        return False, "No aromatic ether bridge linking the two units found"

    # If both conditions are satisfied, we classify the molecule as a bisbenzylisoquinoline alkaloid.
    return True, "Contains at least two benzylisoquinoline units linked by an aromatic ether bridge"