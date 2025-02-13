"""
Classifies: CHEBI:133004 bisbenzylisoquinoline alkaloid
"""
"""
Classifies: Bisbenzylisoquinoline alkaloid
Definition: A benzylisoquinoline alkaloid whose structures are built up of two benzylisoquinoline units 
linked by ether bridges. Various additional bridging patterns (direct C-C linking or methylenedioxy groups)
are common.
"""
from rdkit import Chem

def is_bisbenzylisoquinoline_alkaloid(smiles: str):
    """
    Determines if a molecule is a bisbenzylisoquinoline alkaloid based on its SMILES string.
    Our strategy:
      1. Parse the SMILES string and remove stereochemical annotations to allow for a more permissive
         pattern match.
      2. Use a relaxed SMARTS pattern for a benzylisoquinoline unit. The pattern used here is:
             "c1ccc2CCNC2c1"
         which represents an aromatic ring (c1ccc...) fused to a saturated ring containing nitrogen (CCNC2).
         We then count distinct (non-overlapping) occurrences and require at least two.
      3. Check for the presence of a bridging group between the two units. We look for a generic aromatic
         ether bridge using the SMARTS pattern "[a]-O-[a]" (where [a] stands for any aromatic atom) or a 
         methylenedioxy bridge ("cO[CH2]Oc").
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a bisbenzylisoquinoline alkaloid, False otherwise.
        str: Reason for classification.
    """
    # Parse the SMILES into a molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Remove stereochemistry to allow matching despite chiral annotations.
    Chem.RemoveStereochemistry(mol)
    
    # Define a relaxed SMARTS pattern for a benzylisoquinoline unit.
    # This pattern seeks an aromatic ring fused to a saturated ring containing a nitrogen.
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
    
    # Look for bridging groups.
    # Aromatic ether bridge: an oxygen connecting two aromatic atoms.
    aromatic_ether_smarts = "[a]-O-[a]"
    aromatic_ether = Chem.MolFromSmarts(aromatic_ether_smarts)
    
    # Methylenedioxy bridge: pattern representing a -OCH2O- bridging two aromatic rings.
    methylenedioxy_smarts = "cO[CH2]Oc"
    methylenedioxy = Chem.MolFromSmarts(methylenedioxy_smarts)
    
    has_bridge = mol.HasSubstructMatch(aromatic_ether) or mol.HasSubstructMatch(methylenedioxy)
    if not has_bridge:
        return False, "No bridging group (aromatic ether or methylenedioxy) linking the benzylisoquinoline units found"
    
    # If both conditions are satisfied, classify as bisbenzylisoquinoline alkaloid.
    return True, "Contains at least two benzylisoquinoline units bridged by an aromatic ether or methylenedioxy group"