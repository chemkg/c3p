"""
Classifies: CHEBI:23849 diterpenoid
"""
"""
Classifies: diterpenoid
Definition: Any terpenoid derived from a diterpene. The term includes compounds in which the C20 
skeleton of the parent diterpene has been rearranged or modified by removal of one or more skeletal atoms.
Improved heuristic (empirical):
  1. Molecule must fall in a typical molecular weight range (200–800 Da) and contain carbons.
  2. Extract the Bemis–Murcko scaffold.
  3. Count scaffold carbons, compute the fraction of atoms in the scaffold that are carbons,
     and compute the retention ratio (scaffold carbons / total molecule carbons).
  4. Require that the scaffold carbon count lie between 14 and 25; however if the scaffold is highly fused 
     (4 or more rings) then require at least 16 carbons.
  5. Require that the scaffold is mostly hydrocarbon (carbon fraction ≥ 0.80)
  6. Require that a reasonable fraction of the molecule’s carbons is conserved in the scaffold (retention ratio ≥ 0.50)
  7. Also, if the retention ratio is extremely high (>0.92) then the core is “undecorated” and the molecule may be a false positive.
  8. If any step cannot be satisfied, reject with an explanation.
Note: This heuristic is empirical and some edge cases (both false negatives and false positives) may still occur.
"""

from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold

def is_diterpenoid(smiles: str):
    """
    Determines if a molecule is a diterpenoid based on its SMILES using an improved heuristic.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule is classified as a diterpenoid, False otherwise.
        str: A detailed reason for the classification decision.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check overall molecular weight.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200:
        return False, f"Molecular weight too low ({mol_wt:.1f} Da) to be a diterpenoid"
    if mol_wt > 800:
        return False, f"Molecular weight too high ({mol_wt:.1f} Da) to be a typical diterpenoid"
    
    # Count total carbon atoms in the full molecule.
    total_mol_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if total_mol_carbons == 0:
        return False, "No carbon atoms found in molecule"

    # Extract the Bemis–Murcko scaffold.
    try:
        scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    except Exception as e:
        return False, f"Error computing Murcko scaffold: {str(e)}"
    if scaffold is None or scaffold.GetNumAtoms() == 0:
        return False, "Unable to extract a non-empty Murcko scaffold"
    
    # Count carbon atoms in the scaffold.
    scaffold_carbons = sum(1 for atom in scaffold.GetAtoms() if atom.GetAtomicNum() == 6)
    total_scaffold_atoms = scaffold.GetNumAtoms()
    carbon_fraction = scaffold_carbons / total_scaffold_atoms  # Fraction of scaffold atoms that are carbons.
    
    # Calculate retention ratio: fraction of the molecule's carbons that appear in the scaffold.
    retention_ratio = scaffold_carbons / total_mol_carbons
    
    # Count the number of rings in the scaffold.
    ring_count = rdMolDescriptors.CalcNumRings(scaffold)
    
    # Heuristic checks:
    # 1. Check that the scaffold carbon count is within acceptable bounds.
    #    For all molecules, we allow between 14 and 25 carbons.
    lower_bound = 14
    upper_bound = 25
    # If the scaffold is highly fused (4 or more rings) we demand at least 16 carbons.
    if ring_count >= 4:
        if scaffold_carbons < 16:
            return False, f"Scaffold carbon count is {scaffold_carbons}; with {ring_count} rings, expected at least 16 for diterpenoids"
    
    if not (lower_bound <= scaffold_carbons <= upper_bound):
        return False, f"Scaffold carbon count is {scaffold_carbons}; expected between {lower_bound} and {upper_bound} for diterpenoids"
    
    # 2. Check that the scaffold is mostly hydrocarbon.
    if carbon_fraction < 0.80:
        return False, f"Scaffold carbon fraction too low ({carbon_fraction:.2f}); indicates high heteroatom contribution"
    
    # 3. Ensure a significant fraction of the full molecule’s carbons come from the scaffold.
    if retention_ratio < 0.50:
        return False, f"Low retention ratio ({retention_ratio:.2f}); extensive decorations may obscure the terpene core"
    
    # Optionally, reject molecules with hardly any decoration (retention ratio extremely high).
    if retention_ratio > 0.92:
        return False, f"Retention ratio too high ({retention_ratio:.2f}); core appears nearly undecorated, which is uncommon for diterpenoids"
    
    reason = (f"Scaffold carbon count is {scaffold_carbons} (total scaffold atoms={total_scaffold_atoms}, carbon fraction={carbon_fraction:.2f}), "
              f"{ring_count} rings, retention ratio={retention_ratio:.2f}; molecular weight = {mol_wt:.1f} Da, consistent with a diterpenoid")
    return True, reason

# (Optional) simple test if run as main.
if __name__ == "__main__":
    # Test with an example: Gibberellin A93
    test_smiles = "O1C23C(C(C4OC4C2O)(C1=O)C)C(C56C3CCC(O)(C5)C(C6)=C)C(O)=O"
    result, reason = is_diterpenoid(test_smiles)
    print("Result:", result)
    print("Reason:", reason)