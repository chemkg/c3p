"""
Classifies: CHEBI:133004 bisbenzylisoquinoline alkaloid
"""
"""
Classifies: Bisbenzylisoquinoline alkaloid
Definition: A benzylisoquinoline alkaloid whose structures are built up of two benzylisoquinoline units 
linked by ether bridges. Additional bridging patterns (direct C–C linking or methylenedioxy groups)
are common.
"""
from rdkit import Chem

def is_bisbenzylisoquinoline_alkaloid(smiles: str):
    """
    Determines if a molecule is a bisbenzylisoquinoline alkaloid based on its SMILES string.
    
    Strategy:
      1. Parse and remove stereochemistry to allow flexible matching.
      2. Use a relaxed SMARTS pattern for a benzylisoquinoline unit. We use the pattern:
             "c1ccc2c(c1)CCNC2"
         This pattern represents an aromatic ring fused to a saturated ring containing nitrogen.
         We count distinct (non-overlapping) occurrences and require at least two.
      3. Check for the presence of bridging groups that link the two units.
         We accept either an aromatic ether bridge (pattern "[a]-O-[a]") or a methylenedioxy bridge
         ("cO[CH2]Oc").
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a bisbenzylisoquinoline alkaloid, False otherwise.
        str: Reason for the classification.
    """
    # Parse SMILES string into a molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Remove stereochemical information to be more permissive with the match.
    Chem.RemoveStereochemistry(mol)
    
    # Define a relaxed SMARTS pattern for a benzylisoquinoline unit.
    # This pattern looks for an aromatic ring fused with a saturated ring that contains a nitrogen.
    benzylisoquinoline_smarts = "c1ccc2c(c1)CCNC2"
    benzylisoquinoline_pattern = Chem.MolFromSmarts(benzylisoquinoline_smarts)
    if benzylisoquinoline_pattern is None:
        return False, "Error in benzylisoquinoline SMARTS"
    
    # Get all substructure matches and filter for distinct (non-overlapping) matches.
    iso_matches = mol.GetSubstructMatches(benzylisoquinoline_pattern)
    distinct_units = []
    for match in iso_matches:
        match_set = set(match)
        if not any(match_set & existing for existing in distinct_units):
            distinct_units.append(match_set)
    
    # Require at least 2 distinct benzylisoquinoline units.
    if len(distinct_units) < 2:
        return False, f"Found {len(distinct_units)} benzylisoquinoline unit(s); need at least 2"
    
    # Check for bridging groups linking the two units.
    # Aromatic ether bridge: an oxygen bridging two aromatic atoms.
    aromatic_ether_smarts = "[a]-O-[a]"
    aromatic_ether = Chem.MolFromSmarts(aromatic_ether_smarts)
    
    # Methylenedioxy bridge: pattern representing a -OCH2O- bridge between two aromatic rings.
    methylenedioxy_smarts = "cO[CH2]Oc"
    methylenedioxy = Chem.MolFromSmarts(methylenedioxy_smarts)
    
    has_bridge = mol.HasSubstructMatch(aromatic_ether) or mol.HasSubstructMatch(methylenedioxy)
    if not has_bridge:
        return False, "No bridging group (aromatic ether or methylenedioxy) linking the benzylisoquinoline units found"
    
    # If both conditions are satisfied, we classify the molecule as a bisbenzylisoquinoline alkaloid.
    return True, "Contains at least two benzylisoquinoline units bridged by an appropriate group"