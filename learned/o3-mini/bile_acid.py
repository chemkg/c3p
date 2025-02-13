"""
Classifies: CHEBI:3098 bile acid
"""
#!/usr/bin/env python3
"""
Classifies: Bile acids, defined as hydroxy-5β-cholanic acids (or their glycine/taurine amides or sodium salts).
Note: The classification here is based on heuristic substructure searches:
  - Steroid nucleus: at least four rings (the fused rings of cholesterol/cholan skeleton)
  - Carboxylic acid group: present as “-C(=O)O” (expected at C-24)
  - At least one hydroxy group on a ring (indicative of hydroxylation)
  - Molecular weight in a sensible range for bile acids
This heuristic cannot fully capture the specific stereochemistry (e.g. 5β‑configuration) or the full variation in conjugation.
"""

from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

def is_bile_acid(smiles: str):
    """
    Determines if a molecule is a bile acid based on its SMILES string.
    Bile acids are defined (heuristically) as hydroxy-5β-cholanic acids (or as their glycine/taurine amides or sodium salts)
    that contain a steroid nucleus (a fused set of 4 rings), a carboxylic acid group (typically at C-24) and one or more hydroxy groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is classified as a bile acid, False otherwise.
        str: Reason for classification.
    """
    
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # If the molecule is a salt (multiple fragments), combine them into a single molecule.
    if len(Chem.GetMolFrags(mol)) > 1:
        # Use the largest fragment as the main organic moiety.
        frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
        mol = max(frags, key=lambda m: m.GetNumAtoms())
    
    # Check molecular weight – most bile acids (and their conjugates) are roughly in the 300-1000 Da range.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300 or mol_wt > 1100:
        return False, f"Molecular weight ({mol_wt:.1f} Da) outside expected range for bile acids"
    
    # Check for the presence of at least one carboxylic acid group.
    # Carboxylic acid pattern: carbonyl with a hydroxyl neighbor.
    ca_pattern = Chem.MolFromSmarts("C(=O)[OH]")
    if not mol.HasSubstructMatch(ca_pattern):
        return False, "No carboxylic acid group (-C(=O)O) found"
    
    # Check for at least one hydroxy group on a ring.
    # This pattern matches an –OH group attached to a ring carbon.
    hydroxy_ring_pattern = Chem.MolFromSmarts("[#6;R][OH]")
    hydroxy_matches = mol.GetSubstructMatches(hydroxy_ring_pattern)
    if not hydroxy_matches:
        return False, "No hydroxy group on a ring found"
    
    # Count rings. A steroid nucleus is typically 4 fused rings.
    ring_info = mol.GetRingInfo()
    num_rings = ring_info.NumRings()
    if num_rings < 4:
        return False, f"Found {num_rings} rings; expected at least 4 rings for a steroid nucleus"
    
    # Optionally, one might try to see if the overall topology is similar to that of cholanic acid.
    # This heuristic looks for a fused ring system with 3 six-membered rings plus 1 five-membered ring.
    steroid_pattern = Chem.MolFromSmarts("[$([#6]1[#6][#6][#6][#6]1),"  
                                           "$([#6]1[#6][#6][#6]([#6]1))]")  # very rough...
    # Here we don’t require a match but note that bile acids should have a typical steroid framework.
    steroid_match = mol.HasSubstructMatch(steroid_pattern)
    
    # For bile acids the usual configuration is 5β; however, detecting stereochemistry by SMARTS is challenging.
    # We assume that if the molecule has a steroid nucleus along with the CA group and hydroxy pattern,
    # it is likely a bile acid.
    
    # If we reached here, the molecule meets our heuristic criteria.
    return True, ("Molecule has a steroid nucleus (>=4 rings), carboxylic acid group, and hydroxy groups on rings "
                  "consistent with features expected for a bile acid (note: 5β stereochemistry is not explicitly verified).")

# Example usage (uncomment to test):
# test_smiles = "OC(=O)CC[C@H]([C@@]1([C@@]2([C@]([C@]3([C@@]([C@@]4([C@](CC3)(CC=CC4)[H])C)(CC2)[H])[H])(CC1)[H])C)[H])C"
# result, reason = is_bile_acid(test_smiles)
# print(result, reason)