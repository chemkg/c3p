"""
Classifies: CHEBI:33563 glycolipid
"""
"""
Classifies: glycolipid

Definition: A glycolipid is defined as a molecule that comprises a carbohydrate (sugar)
moiety that is glycosidically linked to a lipid part. Typical examples include 1,2-di-O-acylglycerols
linked to one mono-, di-, or trisaccharide home, but many glycolipids are glycosphingolipids or bacterial glycolipids,
which may contain an amide bond and/or a long aliphatic chain.

Heuristics used in this implementation:
 - Detect a sugar part using both a SMARTS pattern and by scanning the molecule’s ring systems for
   sugar-like features (rings of 5 or 6 atoms containing at least one oxygen and a minimum number of hydroxyl substituents).
 - Look for a lipid part via the presence of ester groups in which the acyl chain (at least 6 carbons) is attached,
   or an amide bond with a similar long chain.
 - Alternatively, if neither ester nor amide is detected, check for a long aliphatic chain (≥8 contiguous carbons).
This is an approximate classification.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def has_sugar(mol):
    """
    Attempt to detect a sugar-like ring.
    
    For each ring of size 5 or 6 that contains at least one oxygen,
    we count the number of ring atoms (typically carbons) that have an -OH substituent.
    If at least two such decorations are found, we consider the ring sugar-like.
    """
    ring_info = mol.GetRingInfo()
    for ring in ring_info.AtomRings():
        if len(ring) in (5, 6):
            ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
            # Check that the ring contains at least one oxygen.
            if any(atom.GetSymbol() == "O" for atom in ring_atoms):
                hydroxyl_count = 0
                for atom in ring_atoms:
                    # Look at neighbors off the ring (or even on the ring) that are –OH groups.
                    for nbr in atom.GetNeighbors():
                        # Avoid counting atoms within the ring.
                        if nbr.GetIdx() in ring:
                            continue
                        if nbr.GetSymbol() == "O":
                            # A hydroxyl oxygen will have at least one hydrogen attached.
                            if any(neighbor.GetSymbol() == "H" for neighbor in nbr.GetNeighbors()):
                                hydroxyl_count += 1
                if hydroxyl_count >= 2:
                    return True
    return False

def is_glycolipid(smiles: str):
    """
    Determines if a molecule is a glycolipid based on its SMILES string.
    
    The improved heuristic:
     - First detects a carbohydrate moiety (via SMARTS and ring analysis).
     - Then checks for the lipid part. The lipid is recognized by either:
         • at least two ester groups where the acyl portion consists of a chain with at least 6 carbons,
         • OR presence of an amide bond linked to a long alkyl chain,
         • OR an isolated long contiguous carbon chain (≥8 C's).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule meets our glycolipid criteria, False otherwise.
        str: Reason for the classification decision.
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Heuristic 1: Detect a sugar moiety.
    # Use a SMARTS pattern for a common pyranose form (flexible, may not match all sugars).
    pyranose_pattern = Chem.MolFromSmarts("OC1C(O)C(O)C(O)C(O)C1")
    sugar_found = mol.HasSubstructMatch(pyranose_pattern)
    # If SMARTS does not match, try a more general ring-based sugar detection.
    if not sugar_found:
        sugar_found = has_sugar(mol)
    
    if not sugar_found:
        return False, "No clear carbohydrate (sugar) moiety detected"
    
    # Heuristic 2: Detect the lipid (acyl) part.
    # Look for ester groups in which the acyl chain has at least 6 contiguous carbons.
    ester_lipid_pattern = Chem.MolFromSmarts("[CX3](=O)OCCCCCC")
    ester_matches = mol.GetSubstructMatches(ester_lipid_pattern)
    n_ester_lipid = len(ester_matches)
    
    # Also look for amide bonds with long chain characteristics.
    amide_lipid_pattern = Chem.MolFromSmarts("NC(=O)CCCCCC")
    amide_matches = mol.GetSubstructMatches(amide_lipid_pattern)
    n_amide_lipid = len(amide_matches)
    
    # Heuristic 3: Check for a long aliphatic chain (a crude indicator).
    long_chain_present = "CCCCCCCC" in smiles  # at least 8 contiguous carbons
    
    # Decision logic:
    if n_ester_lipid >= 2:
        return True, f"Contains sugar moiety and {n_ester_lipid} long-chain ester group(s) (diacyl structure expected)"
    elif n_amide_lipid >= 1 and long_chain_present:
        return True, f"Contains sugar moiety and an amide bond with a long acyl chain"
    elif long_chain_present:
        return True, "Contains sugar moiety and a long aliphatic chain, possible glycolipid variant"
    else:
        return False, "Sugar is present but no clear lipid (diacyl, amide long-chain, or long aliphatic chain) feature detected"

# Example usage:
if __name__ == '__main__':
    # Test example: a simplified glycolipid-like structure (for demonstration)
    test_smiles = "OC1C(O)C(O)C(O)C(O)C1OCC(=O)CCCCCCCCCCCC"
    classification, reason = is_glycolipid(test_smiles)
    print("Is glycolipid?", classification)
    print("Reason:", reason)