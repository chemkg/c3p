"""
Classifies: CHEBI:73011 germacranolide
"""
"""
Classifies: Germacranolide – A sesquiterpene lactone based on a germacrane skeleton.
This updated approach first checks for the presence of a 5‐membered lactone ring 
(γ‐lactone) using a SMARTS pattern. Then the Bemis–Murcko scaffold is computed to isolate 
the core of the molecule. Within that scaffold, we search for a 10‐membered ring that is 
predominantly carbons (heuristic: at least 8 of 10 atoms are carbons) and ensure that the 
scaffold contains roughly 15 carbon atoms (sesquiterpene tropism). These additional criteria 
help resolve some of the formerly misclassified cases.
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
         (with at least 8 carbons) suggesting a germacrane skeleton.
      3. The scaffold should also have a carbon count roughly characteristic of a sesquiterpene (≈15 ± a few).
      4. The molecular weight must be above a minimum threshold for sesquiterpene lactones.
      
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is likely a germacranolide.
        str: Reason for classification or error.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # --- Criterion 1: Check for a 5-membered γ-lactone ring ---
    # The SMARTS pattern below looks for a carbonyl group attached to an oxygen that is in a 5-membered ring.
    lactone_smarts = "[CX3](=O)[OX2r5]"
    lactone_pat = Chem.MolFromSmarts(lactone_smarts)
    if not mol.HasSubstructMatch(lactone_pat):
        return False, "No 5-membered γ-lactone group (ring ester) found."
    
    # --- Criterion 2: Evaluate the core scaffold via Bemis–Murcko ---
    try:
        scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    except Exception as e:
        return False, f"Error obtaining scaffold: {e}"
        
    if scaffold is None:
        return False, "Could not obtain a Bemis–Murcko scaffold."
    
    scaffold_ring_info = scaffold.GetRingInfo()
    scaffold_rings = scaffold_ring_info.AtomRings() if scaffold_ring_info is not None else []
    
    has_10_membered = False
    # Look through rings in the scaffold
    for ring in scaffold_rings:
        if len(ring) == 10:
            carbon_count_ring = sum(1 for idx in ring if scaffold.GetAtomWithIdx(idx).GetAtomicNum() == 6)
            if carbon_count_ring >= 8:
                has_10_membered = True
                break
    if not has_10_membered:
        return False, "No suitable 10-membered carbocycle (germacrane skeleton) found in scaffold."
    
    # --- Criterion 3: Check the carbon count in the scaffold (expect about 15 carbons for sesquiterpenes) ---
    scaffold_atoms = scaffold.GetAtoms()
    scaffold_carbon_count = sum(1 for atom in scaffold_atoms if atom.GetAtomicNum() == 6)
    # Allow a tolerance (e.g., 12 to 18 carbons) because substituents might add extra atoms.
    if not (12 <= scaffold_carbon_count <= 18):
        return False, f"Scaffold carbon count ({scaffold_carbon_count}) is not in the expected range for a sesquiterpene."
    
    # --- Criterion 4: Check molecular weight to be above a lower limit, e.g. 200 Da ---
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200:
        return False, "Molecular weight is too low to be a sesquiterpene lactone."
    
    return True, "Contains a 10-membered germacrane skeleton (in the core scaffold), a 5-membered γ-lactone, and appropriate carbon count typical of germacranolides."

# Example use:
if __name__ == "__main__":
    example_smiles = "C\\C=C(/C)C(=O)O[C@H]1C\\C(C)=C\\C(=O)\\C=C(C)/[C@H](O)[C@H]2OC(=O)C(=C)[C@H]12"  # Molephantinin
    result, reason = is_germacranolide(example_smiles)
    print(result, reason)