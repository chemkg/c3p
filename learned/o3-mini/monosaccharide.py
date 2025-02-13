"""
Classifies: CHEBI:35381 monosaccharide
"""
"""
Classifies: Monosaccharide
A monosaccharide is defined here as a single (non‐glycosidically connected) polyhydroxy aldehyde or ketone 
(with potential to exist as a cyclic hemiacetal), having three or more carbon atoms. 
This heuristic checks for a single fragment, at least three carbons, the presence of a carbonyl (for open‐chain forms)
or a furanose/pyranose–like ring (for cyclic forms) and sufficient hydroxyl groups.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_monosaccharide(smiles: str):
    """
    Determines if a molecule is a monosaccharide based on its SMILES string.
    The compound must be a single unit (no glycosidic connections) with at least three carbons,
    multiple hydroxyl groups, and either an explicit carbonyl (aldehyde or ketone) for open-chain forms
    or a sugar–like ring (furanose or pyranose) for cyclic forms.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is classified as a monosaccharide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check that molecule is a single fragment (i.e. one connected unit)
    frags = Chem.GetMolFrags(mol, asMols=True)
    if len(frags) > 1:
        return False, "Multiple fragments detected – likely a glycosidic conjugate"
    
    # Count carbon atoms; a monosaccharide should have at least 3 carbons.
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 3:
        return False, f"Too few carbons ({carbon_count} found, need at least 3)"
    
    # Add explicit hydrogens to help find –OH groups in substructure search
    molH = Chem.AddHs(mol)
    
    # Count hydroxyl (–OH) groups.
    # The pattern [OX2H] finds oxygen with two connections that has an attached hydrogen.
    OH_pattern = Chem.MolFromSmarts("[OX2H]")
    hydroxyl_matches = molH.GetSubstructMatches(OH_pattern)
    hydroxyl_count = len(hydroxyl_matches)
    if hydroxyl_count < 2:
        return False, f"Too few hydroxyl groups ({hydroxyl_count} found), expected a polyhydroxy compound"
    
    # Look for an explicit carbonyl group.
    # This SMARTS will match a C=O group (which can be part of an aldehyde or ketone).
    carbonyl_pattern = Chem.MolFromSmarts("[CX3](=O)")
    has_carbonyl = mol.HasSubstructMatch(carbonyl_pattern)
    
    # Look for a sugar-like ring: furanose (5-membered) or pyranose (6-membered)
    has_sugar_ring = False
    ring_info = mol.GetRingInfo()
    for ring in ring_info.AtomRings():
        # Consider only rings of size 5 or 6.
        if len(ring) in (5, 6):
            # Check how many oxygen atoms are in the ring.
            oxygens_in_ring = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8)
            # In sugars one of the atoms is the ring-oxygen.
            if oxygens_in_ring == 1:
                has_sugar_ring = True
                break
    
    # Evaluate:
    # A monosaccharide parent unit will either have a free carbonyl (open-chain) or exist as a cyclic sugar.
    if not has_carbonyl and not has_sugar_ring:
        return False, "No carbonyl or typical sugar ring detected; not a monosaccharide motif"
    
    # If we found a carbonyl, make sure it is not part of an unrelated function.
    # (Advanced checks could be included here but are omitted for simplicity.)
    
    # Passed our heuristic tests.
    reason = ("Molecule is a single unit with {} carbons, {} hydroxyl groups, and {} detected."
              .format(carbon_count, hydroxyl_count,
                      "an explicit carbonyl" if has_carbonyl else "a sugar-like ring"))
    return True, reason

# For testing, you can uncomment examples below:
# examples = [
#     "O1[C@@]([C@H](O)[C@@H](O)[C@@H]1O)([C@H](O)CO)[H]",  # beta-D-glucofuranose
#     "[H]C(=O)[C@H](O)[C@@H](O)[C@@H](O)CO",              # aldehydo-L-arabinose
#     "OC[C@@H]1O[C@@H](O)[C@@H](O)[C@H](O)[C@@H]1O",       # alpha-L-galactose
# ]
# for ex in examples:
#     is_sugar, msg = is_monosaccharide(ex)
#     print(ex, is_sugar, msg)