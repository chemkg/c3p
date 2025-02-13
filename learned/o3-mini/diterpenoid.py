"""
Classifies: CHEBI:23849 diterpenoid
"""
"""
Classifies: diterpenoid
Definition: Any terpenoid derived from a diterpene. The term includes compounds in which the C20 
skeleton of the parent diterpene has been rearranged or modified by the removal of one or more skeletal atoms.
Revised heuristic:
  1. First, check that the overall molecular weight is in a range typical for diterpenoids (e.g., 200–800 Da).
  2. Extract the Bemis–Murcko scaffold.
  3. Compute the carbon count of the scaffold and the fraction of heavy atoms in the scaffold that are carbon.
  4. Also compute the ratio of scaffold carbons to the molecule’s total carbon count (“retention ratio”).
  5. Accept if:
       • The scaffold carbon count is within an acceptable range (we use 14–25; allowing for loss from C20),
         but if the scaffold is very fused (ring count > 3) then require that the count is close to 20 (e.g. 18–22).
       • The scaffold’s carbon fraction is high (>= 0.80, meaning the core is mostly hydrocarbon).
       • The retention ratio (scaffold carbon count / total molecular carbon count) is at least 0.5,
         meaning that a good part of the molecule comes from the original terpene core.
  6. Otherwise, reject.
Note:
  This heuristic is empirical and may mis‐classify some edge cases.
"""

from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold

def is_diterpenoid(smiles: str):
    """
    Determines if a molecule is a diterpenoid based on its SMILES string using an improved heuristic.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is classified as a diterpenoid, False otherwise.
        str: Reason for classification.
    """
    
    # Parse SMILES into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check overall molecular weight (diterpenoids typically are not very small or huge).
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200:
        return False, f"Molecular weight too low ({mol_wt:.1f} Da) to be a diterpenoid"
    if mol_wt > 800:
        return False, f"Molecular weight too high ({mol_wt:.1f} Da) to be a typical diterpenoid"
    
    # Count overall carbon atoms in the full molecule.
    total_mol_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if total_mol_carbons == 0:
        return False, "No carbon atoms found in molecule"
    
    # Try to extract the Bemis–Murcko scaffold.
    try:
        scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    except Exception as e:
        return False, f"Error computing Murcko scaffold: {str(e)}"
    if scaffold is None or scaffold.GetNumAtoms() == 0:
        return False, "Unable to extract a non-empty Murcko scaffold"
    
    # Count the number of carbon atoms in the scaffold.
    scaffold_carbons = sum(1 for atom in scaffold.GetAtoms() if atom.GetAtomicNum() == 6)
    total_scaffold_atoms = scaffold.GetNumAtoms()
    carbon_fraction = scaffold_carbons / total_scaffold_atoms  # fraction of scaffold atoms that are carbon
    
    # Calculate retention ratio: how much of the full molecule’s carbons appear in the scaffold.
    retention_ratio = scaffold_carbons / total_mol_carbons
    
    # Count number of rings in the scaffold.
    ring_count = rdMolDescriptors.CalcNumRings(scaffold)
    
    # Now apply our criteria.
    # Set acceptable scaffold carbon count range. Note: we allow a bit of variation from C20.
    if ring_count > 3:
        # If highly fused, require the carbon count be closer to C20.
        lower_bound = 18
        upper_bound = 22
    else:
        lower_bound = 14
        upper_bound = 25
        
    # Check scaffold carbon count.
    if not (lower_bound <= scaffold_carbons <= upper_bound):
        return (False, f"Scaffold carbon count is {scaffold_carbons}; expected between {lower_bound} and {upper_bound} for diterpenoids")
    
    # Check that the scaffold is mostly hydrocarbon.
    if carbon_fraction < 0.80:
        return (False, f"Scaffold carbon fraction too low ({carbon_fraction:.2f}); indicates heavy heteroatom content")
    
    # Ensure that a significant fraction of the molecule comes from the original terpene core.
    if retention_ratio < 0.50:
        return (False, f"Low retention ratio ({retention_ratio:.2f}); core structure may be obscured by decorations")
    
    # If we reach here, then the heuristic considers it consistent with a diterpenoid.
    reason = (f"Scaffold carbon count is {scaffold_carbons} (total atoms={total_scaffold_atoms}, carbon fraction={carbon_fraction:.2f}), "
              f"{ring_count} rings, retention ratio={retention_ratio:.2f}; molecular weight = {mol_wt:.1f} Da, "
              "consistent with a diterpenoid")
    return True, reason

# (Optional) simple test when running module stand-alone.
if __name__ == "__main__":
    # Test with one example: Gibberellin A93
    test_smiles = "O1C23C(C(C4OC4C2O)(C1=O)C)C(C56C3CCC(O)(C5)C(C6)=C)C(O)=O"
    result, reason = is_diterpenoid(test_smiles)
    print("Result:", result)
    print("Reason:", reason)