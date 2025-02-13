"""
Classifies: CHEBI:73011 germacranolide
"""
"""
Classifies: Germacranolide – A sesquiterpene lactone based on a germacrane skeleton.
Criteria used here (heuristic):
  1. The molecule must contain a 5‐membered lactone ring (γ‐lactone). We search for
     a ring ester using the SMARTS "[CX3](=O)[OX2r5]".
  2. The molecule should have a 10‑membered carbocycle (the germacrane skeleton).
     We inspect the ring information and require a ring of size 10 where the vast majority
     of atoms are carbons.
  3. As a loose check, we also verify that the molecular weight is in a reasonable range.
  
Note: This classification is heuristic and may not catch all edge cases.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_germacranolide(smiles: str):
    """
    Determines if a molecule is a germacranolide based on its SMILES string.
    A germacranolide is a sesquiterpene lactone based on a germacrane skeleton.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule features a 10-membered germacrane skeleton and
              a γ-lactone (5-membered lactone ring), typical for germacranolides.
        str: Reason for the classification.
    """
    # Parse SMILES and check if valid
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # --- Check for a 5-membered lactone ring (γ-lactone) ---
    # The SMARTS below searches for a carbonyl carbon (C(=O)) attached to an oxygen that is
    # part of a ring (r5 indicates a ring atom in a 5-membered ring).
    lactone_smarts = "[CX3](=O)[OX2r5]"
    lactone_pat = Chem.MolFromSmarts(lactone_smarts)
    if not mol.HasSubstructMatch(lactone_pat):
        return False, "No 5-membered γ-lactone group (ring ester) found."
    
    # --- Check for a 10-membered carbocycle (germacrane skeleton) ---
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    has_10_membered = False
    for ring in rings:
        if len(ring) == 10:
            # Count the number of carbon atoms in this ring
            carbon_count = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
            # Require that most atoms in the ring are carbons (heuristic: at least 8 out of 10)
            if carbon_count >= 8:
                has_10_membered = True
                break
    if not has_10_membered:
        return False, "No suitable 10-membered carbocycle (germacrane skeleton) found."
    
    # --- Optional: Check if the molecular weight is in a plausible range ---
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200:
        return False, "Molecular weight is too low to be a sesquiterpene lactone."
    
    # If all criteria are met, classify as a germacranolide.
    return True, "Contains a 10-membered germacrane skeleton and a 5-membered γ-lactone group typical of germacranolides."