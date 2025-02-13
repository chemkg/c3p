"""
Classifies: CHEBI:35341 steroid
"""
#!/usr/bin/env python
"""
Classifies: Steroid
Definition:
  Any of naturally occurring compounds and synthetic analogues, based on the cyclopenta[a]phenanthrene carbon skeleton,
  partially or completely hydrogenated; there are usually methyl groups at C-10 and C-13, and often an alkyl group at C-17.
  By extension, one or more bond scissions, ring expansions and/or ring contractions of the skeleton may have occurred.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_steroid(smiles: str):
    """
    Determines if a molecule is a steroid based on its SMILES string.
    The classification uses two strategies:
      1. A primary substructure search for a canonical cyclopenta[a]phenanthrene nucleus.
      2. If the primary match fails, a secondary fuzzy check is performed by analyzing the fused ring system --
         requiring at least three six-membered rings and at least one five-membered ring.
         
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a steroid, False otherwise.
        str: Reason for the classification decision.
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Primary check: Look for the canonical steroid nucleus.
    # The SMARTS "C1CCC2C3CCC4C1C3CCC2C4" represents the core of the cyclopenta[a]phenanthrene skeleton.
    steroid_pattern = Chem.MolFromSmarts("C1CCC2C3CCC4C1C3CCC2C4")
    if mol.HasSubstructMatch(steroid_pattern):
        return True, "Molecule contains the canonical steroid nucleus (cyclopenta[a]phenanthrene skeleton)."
    
    # Secondary check: Analyze the ring system in the molecule.
    # Retrieve ring information from the molecule.
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()  # Each ring is a tuple of atom indices.
    
    # Count rings by size.
    six_membered = sum(1 for ring in rings if len(ring) == 6)
    five_membered = sum(1 for ring in rings if len(ring) == 5)
    
    # For a steroid nucleus one expects three six-membered rings and one five-membered ring.
    if six_membered >= 3 and five_membered >= 1:
        return True, ("Molecule has a fused ring system with at least three six-membered and one five-membered ring, "
                      "suggestive of a steroid nucleus (possibly modified).")
    
    return False, "Molecule does not appear to contain a steroid-like fused ring system."

# Example test cases (you can uncomment for testing):
# smiles_examples = [
#     "O=C1C[C@]2([C@@]([C@@]3(C([C@]4([C@@]([C@](CC4)([C@@H](CC/C(/C(C)C)=C/C)C)[H])(CC3)C)[H])=CC2)[H])(CC1)C",  # Avenastenone (steroid)
#     "O=C1NC(=NC2=C1N(C=N2)C[C@]3(O)[C@@]4(OC=5C6=C(O)C=7C(=O)C(O)CC(C7C(=C6C=C(C5[C@H]8[C@@H]4O[C@@]3(O8)C(OC)OC)C)OC)O[C@@H]9OC([C@@H](OC(=O)C)C(C9)(O)C)C)O[C@@H]%10OC([C@](O)(C(=O)C)C(C%10)O)C)N",  # Gutingimycin (non-steroid)
# ]
#
# for s in smiles_examples:
#     result, reason = is_steroid(s)
#     print("Result:", result, "\nReason:", reason, "\n")