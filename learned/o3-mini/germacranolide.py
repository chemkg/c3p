"""
Classifies: CHEBI:73011 germacranolide
"""
"""
Classifies: Germacranolide – A sesquiterpene lactone based on a germacrane skeleton.
This updated approach first checks for the presence of a 5‐membered lactone ring 
(γ‐lactone) using a SMARTS pattern. Then the Bemis–Murcko scaffold is computed to 
isolate the core of the molecule. Within that scaffold (or in the parent molecule if
the scaffold fails to exhibit a clear 10‐membered ring), we search for a 10‐membered ring 
that is predominantly carbons (heuristic: at least 8 out of 10 atoms) and ensure that the 
scaffold exhibits a carbon count typical for sesquiterpenes (roughly 10–20 carbons).
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold

def is_germacranolide(smiles: str):
    """
    Determines if a molecule is a germacranolide based on its SMILES string.
    A germacranolide is a sesquiterpene lactone based on a germacrane skeleton.
    
    The classification now uses the following heuristic criteria:
      1. The molecule must contain a 5‐membered lactone ring (γ‐lactone), detected using the SMARTS "[CX3](=O)[OX2r5]".
      2. The molecule’s core scaffold (Bemis–Murcko scaffold) should contain a 10‐membered ring 
         (with at least 8 carbons) suggesting a germacrane skeleton. If no such ring is found,
         then the original molecule is examined as a fallback.
      3. If a scaffold is available, its carbon count should be in the range of roughly 10 to 20 carbons,
         typical for sesquiterpenes.
      4. The molecular weight must be above 200 Da.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule is likely a germacranolide, False otherwise.
        str: Reason for classification.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # Criterion 1: Check for a 5-membered γ-lactone ring using SMARTS pattern.
    lactone_smarts = "[CX3](=O)[OX2r5]"
    lactone_pat = Chem.MolFromSmarts(lactone_smarts)
    if not mol.HasSubstructMatch(lactone_pat):
        return False, "No 5-membered γ-lactone group (ring ester) found."
    
    # Criterion 2: Get the Bemis–Murcko scaffold.
    scaffold = None
    try:
        scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    except Exception as e:
        # If scaffold cannot be calculated, we continue using the molecule itself.
        scaffold = None
    
    has_10_membered = False
    
    # First, try to find a 10-membered ring in the scaffold (if available).
    if scaffold is not None:
        ring_info = scaffold.GetRingInfo()
        rings = ring_info.AtomRings() if ring_info is not None else []
        for ring in rings:
            if len(ring) == 10:
                # Count carbons in this ring.
                carbon_count = sum(1 for idx in ring 
                                   if scaffold.GetAtomWithIdx(idx).GetAtomicNum() == 6)
                if carbon_count >= 8:
                    has_10_membered = True
                    break
    
    # Fallback: if no suitable ring was found in the scaffold, search the full molecule.
    if not has_10_membered:
        ring_info = mol.GetRingInfo()
        rings = ring_info.AtomRings() if ring_info is not None else []
        for ring in rings:
            if len(ring) == 10:
                carbon_count = sum(1 for idx in ring
                                   if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
                if carbon_count >= 8:
                    has_10_membered = True
                    break
    
    if not has_10_membered:
        return False, "No suitable 10-membered carbocycle (germacrane skeleton) found."
    
    # Criterion 3: If a scaffold is available, check if its carbon count is roughly in the expected range.
    if scaffold is not None:
        scaffold_carbon_count = sum(1 for atom in scaffold.GetAtoms() if atom.GetAtomicNum() == 6)
        if not (10 <= scaffold_carbon_count <= 20):
            return False, f"Scaffold carbon count ({scaffold_carbon_count}) is not in the expected range for sesquiterpenes."
    
    # Criterion 4: Check if the molecular weight exceeds the minimal threshold.
    mw = rdMolDescriptors.CalcExactMolWt(mol)
    if mw < 200:
        return False, "Molecular weight is too low to be a sesquiterpene lactone."
    
    return True, ("Contains a 10-membered germacrane skeleton (in the core scaffold or molecule), "
                  "a 5-membered γ-lactone, and an appropriate carbon count typical of germacranolides.")

# Example usage:
if __name__ == "__main__":
    # Example with Molephantinin SMILES:
    example_smiles = "C\\C=C(/C)C(=O)O[C@H]1C\\C(C)=C\\C(=O)\\C=C(C)/[C@H](O)[C@H]2OC(=O)C(=C)[C@H]12"
    result, reason = is_germacranolide(example_smiles)
    print(result, reason)