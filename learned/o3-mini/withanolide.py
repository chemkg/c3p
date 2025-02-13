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
    A withanolide is defined as a steroid lactone having a C28 steroid nucleus and a lactone ring
    formed from a modified side chain, including its substituted derivatives.
    
    The criteria used here are heuristic:
      - The molecule must be valid.
      - It should contain a steroid nucleus. Here we use a generic SMARTS for the cyclopentanoperhydrophenanthrene core.
      - It should contain at least one lactone ring (cyclic ester).
      - The molecule should contain at least 28 carbons.
      - It is expected to have at least 5 rings overall (4 from the steroid nucleus and 1 from the lactone).
      
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule is classified as withanolide, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for a steroid nucleus.
    # The steroid nucleus (cyclopentanoperhydrophenanthrene) typically has 3 six-membered rings and 1 five-membered ring fused.
    # This is a simplistic SMARTS pattern representing a fused tetracyclic system:
    steroid_pattern = Chem.MolFromSmarts("C1CC2CCC3C(C2)C1CCC3")
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid nucleus found"
        
    # Count carbon atoms.
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 28:
        return False, f"Only {c_count} carbon atoms found, fewer than expected for a C28 steroid nucleus"
    
    # Check overall ring count.
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 5:
        return False, f"Only {ring_info.NumRings()} rings found, expected at least 5 rings (steroid nucleus + lactone ring)"
        
    # Look for a lactone (cyclic ester) ring.
    # We use a SMARTS pattern for a carbonyl (C(=O)) bound to an oxygen which is part of a ring.
    lactone_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2r]")
    if not mol.HasSubstructMatch(lactone_pattern):
        return False, "No lactone ring (cyclic ester) found"
    
    # If all criteria are met, the molecule is classified as withanolide.
    return True, "Contains a steroid nucleus with at least 28 carbons and a lactone ring, consistent with withanolide structure"

# Example usage (uncomment to test):
# smiles_example = "C[C@@H]([C@H]1CC[C@H]2[C@@H]3C[C@H]4O[C@]44[C@@H](O)C=CC(=O)[C@]4(C)[C@H]3CC[C@]12C)[C@H]1CC(C)=C(CO)C(=O)O1"  # (withaferin A-like)
# result, reason = is_withanolide(smiles_example)
# print(result, reason)