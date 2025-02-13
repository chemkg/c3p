"""
Classifies: CHEBI:23849 diterpenoid
"""
"""
Classifies: diterpenoid
Definition: Any terpenoid derived from a diterpene. The term includes compounds in which the C20 
skeleton of the parent diterpene has been rearranged or modified by the removal of one or more skeletal atoms
(generally methyl groups).

Heuristic: We extract the Bemis–Murcko scaffold and count the number of carbon atoms.
If the scaffold contains between 16 and 23 carbon atoms (a range that allows for small losses from C20),
and the molecule is not extremely light, we label it as a diterpenoid.
"""

from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold

def is_diterpenoid(smiles: str):
    """
    Determines if a molecule is a diterpenoid based on its SMILES string.
    Using a heuristic based on the Bemis–Murcko scaffold: diterpenoids are derived from a C20 terpene
    (possibly losing one or more carbon atoms by rearrangement).
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is classified as a diterpenoid, False otherwise.
        str: Reason for classification.
    """
    # Parse SMILES into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Calculate molecular weight as a preliminary check (very low molecular weight unlikely to be a diterpenoid)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200:
        return False, f"Molecular weight too low ({mol_wt:.1f} Da) to be a diterpenoid"
    
    # Get the Bemis–Murcko scaffold as the core of the molecule.
    try:
        scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    except Exception as e:
        return False, f"Error computing Murcko scaffold: {str(e)}"
    
    if scaffold is None:
        return False, "Unable to extract Murcko scaffold"

    # Count the number of carbon atoms in the scaffold
    c_count = sum(1 for atom in scaffold.GetAtoms() if atom.GetAtomicNum() == 6)
    
    # Use a heuristic range: since the diterpene parent is C20, allowance is made for small losses.
    # We use an approximate range of 16 to 23 carbons in the core scaffold.
    if 16 <= c_count <= 23:
        return True, f"Scaffold carbon count is {c_count}, which is consistent with a diterpenoid"
    else:
        return False, f"Scaffold carbon count is {c_count}; not in the expected range for diterpenoids (16-23 C's)"
        
# (Optional) For testing purposes, you could call the function with one of the provided SMILES:
if __name__ == "__main__":
    test_smiles = "O1C23C(C(C4OC4C2O)(C1=O)C)C(C56C3CCC(O)(C5)C(C6)=C)C(O)=O"  # gibberellin A93 example
    result, reason = is_diterpenoid(test_smiles)
    print("Result:", result)
    print("Reason:", reason)