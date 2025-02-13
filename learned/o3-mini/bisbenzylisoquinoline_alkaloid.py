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
    The classification uses two main criteria:
      1. The molecule must contain at least two distinct benzylisoquinoline units.
         We search for the bicyclic fragment representing a tetrahydroisoquinoline unit with an attached aromatic ring.
         A more permissive SMARTS pattern "c1ccc2CCNC2c1" is used to allow for additional substituents and loss
         of explicit stereochemistry.
      2. The two units must be connected by a bridging group. We check for:
         - An aromatic ether bridge (SMARTS: "[c;r]-O-[c;r]"), or
         - A methylenedioxy bridge (SMARTS: "cO[CH2]Oc").
    
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
    
    # Define a more permissive SMARTS pattern for a benzylisoquinoline unit.
    # This pattern matches an aromatic ring fused to a saturated ring containing nitrogen.
    benzylisoquinoline_smarts = "c1ccc2CCNC2c1"
    benzylisoquinoline_pattern = Chem.MolFromSmarts(benzylisoquinoline_smarts)
    iso_matches = mol.GetSubstructMatches(benzylisoquinoline_pattern)
    
    # Filter matches to obtain distinct (non-overlapping) units.
    distinct_units = []
    for match in iso_matches:
        match_set = set(match)
        if not any(match_set & existing for existing in distinct_units):
            distinct_units.append(match_set)
    
    if len(distinct_units) < 2:
        return False, f"Found {len(distinct_units)} benzylisoquinoline unit(s); need at least 2"
    
    # Look for bridging patterns.
    # Aromatic ether bridge: an oxygen connecting two aromatic (ring) carbons.
    ether_bridge_smarts = "[c;r]-O-[c;r]"
    ether_bridge = Chem.MolFromSmarts(ether_bridge_smarts)
    
    # Methylenedioxy bridge: pattern representing a -OCH2O- group bridging two aromatic rings.
    methylenedioxy_smarts = "cO[CH2]Oc"
    methylenedioxy = Chem.MolFromSmarts(methylenedioxy_smarts)
    
    has_bridge = mol.HasSubstructMatch(ether_bridge) or mol.HasSubstructMatch(methylenedioxy)
    if not has_bridge:
        return False, "No bridging group (aromatic ether or methylenedioxy) linking the benzylisoquinoline units found"
    
    # If both conditions are satisfied, the molecule is classified as a bisbenzylisoquinoline alkaloid.
    return True, "Contains at least two benzylisoquinoline units bridged by an aromatic ether or methylenedioxy group"