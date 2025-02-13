"""
Classifies: CHEBI:23849 diterpenoid
"""
"""
Classifies: diterpenoid
Definition: Any terpenoid derived from a diterpene. The term includes compounds in which the C20 
skeleton of the parent diterpene has been rearranged or modified by the removal of one or more skeletal atoms.
Heuristic:
  1. Extract the Bemis–Murcko scaffold and count the number of carbon atoms.
  2. Accept if the scaffold carbon count is between 15 and 23 (allowing for some loss from C20),
     the fraction of carbons in the scaffold is high (>=0.75), 
     and the scaffold does not include too many rings (<=3 rings).
  3. In addition, reject molecules with very low molecular weight.
"""

from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold

def is_diterpenoid(smiles: str):
    """
    Determines if a molecule is a diterpenoid based on its SMILES string.
    Using a heuristic based on the Bemis–Murcko scaffold (carbon count, heteroatom content,
    and ring count). Diterpenoids are derived from a nominal C20 terpene so we allow a modest loss
    in carbons (accept scaffold carbon counts between 15 and 23) and require that the core is mostly
    hydrocarbon (carbon fraction >= 0.75) and not heavily aromatic or overly fused (<=3 rings in the scaffold).
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is classified as a diterpenoid, False otherwise.
        str: Reason for classification.
    """
    # Parse the SMILES string into a molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Preliminary check: very low molecular weight compounds are unlikely to be diterpenoids.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200:
        return False, f"Molecular weight too low ({mol_wt:.1f} Da) to be a diterpenoid"
    
    # Extract the Bemis–Murcko scaffold for the molecule.
    try:
        scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    except Exception as e:
        return False, f"Error computing Murcko scaffold: {str(e)}"
    if scaffold is None:
        return False, "Unable to extract Murcko scaffold"
    
    # Count how many heavy atoms in the scaffold are carbon (atomic number 6)
    c_count = sum(1 for atom in scaffold.GetAtoms() if atom.GetAtomicNum() == 6)
    total_atoms = scaffold.GetNumAtoms()
    if total_atoms == 0:
        return False, "Empty scaffold"
    carbon_ratio = c_count / total_atoms  # fraction of heavy atoms that are carbon

    # Count the number of rings in the scaffold.
    ring_count = rdMolDescriptors.CalcNumRings(scaffold)
    
    # Check the heuristic criteria.
    # (1) The basic carbon count: for a diterpene-derived structure we expect ~C20, 
    #     with modifications giving a core of 15 to 23 carbon atoms.
    if not (15 <= c_count <= 23):
        return False, f"Scaffold carbon count is {c_count}; not in the expected range for diterpenoids (15-23 C's)"
    
    # (2) The scaffold should be mostly hydrocarbon.
    if carbon_ratio < 0.75:
        return False, f"Scaffold carbon fraction too low ({carbon_ratio:.2f}); indicates heavy heteroatom content"
    
    # (3) A fused ring system that is very dense (e.g., four or more rings) is more typical of steroids.
    if ring_count > 3:
        return False, f"Scaffold has {ring_count} rings; too many rings for a typical diterpenoid scaffold"
    
    # Passed all checks: classify as diterpenoid.
    return True, f"Scaffold carbon count is {c_count} (ratio={carbon_ratio:.2f}) with {ring_count} rings, consistent with a diterpenoid"

# (Optional) Simple test when running this module stand-alone.
if __name__ == "__main__":
    # Test with one example: Gibberellin A93
    test_smiles = "O1C23C(C(C4OC4C2O)(C1=O)C)C(C56C3CCC(O)(C5)C(C6)=C)C(O)=O"
    result, reason = is_diterpenoid(test_smiles)
    print("Result:", result)
    print("Reason:", reason)