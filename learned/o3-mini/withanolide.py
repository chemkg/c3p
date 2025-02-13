"""
Classifies: CHEBI:74716 withanolide
"""
"""
Classifies: Withanolide – any steroid lactone that is a C28 steroid with a modified side chain forming a lactone ring and its substituted derivatives.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_withanolide(smiles: str):
    """
    Determines if a molecule is a withanolide based on its SMILES string.
    A withanolide is defined as a steroid lactone having a C28 steroid nucleus (or modified derivative)
    and a lactone ring formed from a modified side chain.
    
    The criteria used here (heuristic) are as follows:
      - The SMILES should be valid.
      - The molecule must have at least 28 carbon atoms.
      - The molecule should have at least 4 rings overall.
      - Among its rings, at least one should have 5 atoms and at least three should have 6 atoms.
        (This is a proxy for the fused steroid nucleus – typically the steroid core is composed of three six‐membered rings and one five‐membered ring fused together.)
      - The molecule must contain at least one lactone ring (a cyclic ester), identified via a SMARTS match.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule is classified as withanolide, False otherwise.
        str: Reason for the classification decision.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Count carbon atoms.
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 28:
        return False, f"Only {c_count} carbon atoms detected, fewer than expected for a C28 steroid nucleus"
    
    # Get the ring information from the molecule
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()  # tuple of tuples containing atom indices in each ring
    n_rings = len(rings)
    
    if n_rings < 4:
        return False, f"Only {n_rings} rings found, expected at least 4 rings (steroid nucleus + lactone ring)"
    
    # Heuristic: in a steroid, we expect one five-membered ring and three six-membered rings
    ring_sizes = [len(r) for r in rings]
    n_five = sum(1 for size in ring_sizes if size == 5)
    n_six = sum(1 for size in ring_sizes if size == 6)
    if n_five < 1 or n_six < 3:
        return False, (f"Rings detected: sizes {ring_sizes}. Expected at least one 5-membered ring and three 6-membered rings "
                       "characteristic of a steroid nucleus")
    
    # Look for a lactone (cyclic ester) ring.
    # This SMARTS pattern seeks a carbonyl group (C(=O)) attached to an oxygen that is part of a ring.
    lactone_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2r]")
    if not mol.HasSubstructMatch(lactone_pattern):
        return False, "No lactone ring (cyclic ester) found"
    
    # If all criteria are met, then classify as withanolide.
    return True, ("Molecule has at least 28 carbons, a fused ring system consistent with a modified steroid nucleus "
                  "and a lactone ring, consistent with a withanolide structure.")

# Example usage (uncomment to test):
# smiles_example = "C[C@@H]([C@H]1CC[C@H]2[C@@H]3C[C@H]4O[C@]44[C@@H](O)C=CC(=O)[C@]4(C)[C@H]3CC[C@]12C)[C@H]1CC(C)=C(CO)C(=O)O1"  # (withaferin A-like)
# result, reason = is_withanolide(smiles_example)
# print(result, reason)