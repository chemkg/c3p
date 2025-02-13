"""
Classifies: CHEBI:39434 limonoid
"""
"""
Classifies: Limonoid
Definition:
    Any triterpenoid that is highly oxygenated and has a prototypical structure either containing or derived from a
    precursor with a 4,4,8-trimethyl-17-furanylsteroid skeleton.
    (e.g. limonin, azadirachtin, and other structurally related compounds)
    
This implementation uses heuristic criteria:
    - The molecule must be sufficiently large (MW > 300 and enough carbon atoms, here we require at least 20)
    - It must be “highly oxygenated” (at least 4 oxygen atoms)
    - It must be polycyclic (at least 4 rings)
    - It must contain a furan ring (SMARTS "c1ccoc1")
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_limonoid(smiles: str):
    """
    Determines if a molecule is a limonoid based on its SMILES string.
    
    Strategy (heuristic):
      - Check that the SMILES can be parsed.
      - Verify molecular weight > 300 Da and at least 20 carbon atoms (suggestive of a triterpenoid framework).
      - Ensure the structure is highly oxygenated (a minimum number of oxygen atoms).
      - Verify the molecule has a polycyclic core (at least 4 rings).
      - Ensure that a furan ring is present (using a SMARTS search for "c1ccoc1").
      
    Args:
       smiles (str): SMILES representation of the molecule.
       
    Returns:
       bool: True if the molecule is classified as a limonoid, otherwise False.
       str: Description of the reasoning.
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check molecular weight (limonoids are derived from triterpenoids and typically are >300 Da)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:
        return False, f"Molecular weight is too low ({mol_wt:.1f} Da); expected >300 Da for a triterpenoid derivative"
    
    # Check number of carbon atoms: triterpenoid-derived molecules must have a high carbon count
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 20:
        return False, f"Too few carbon atoms ({c_count}); expected at least 20 for a triterpenoid"
    
    # Check for high oxygenation: count oxygen atoms (heuristically require at least 4)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 4:
        return False, f"Not enough oxygen atoms ({o_count}); limonoids are highly oxygenated"
    
    # Check for polycyclic structure: we expect at least 4 rings
    ring_info = mol.GetRingInfo()
    n_rings = ring_info.NumRings()
    if n_rings < 4:
        return False, f"Not enough cyclic rings ({n_rings}); expected a polycyclic structure"
    
    # Check for the characteristic furan ring using a SMARTS pattern
    furan_pattern = Chem.MolFromSmarts("c1ccoc1")
    if not mol.HasSubstructMatch(furan_pattern):
        return False, "Missing furan ring which is typical in limonoids"
    
    return True, "Molecule appears to be a limonoid (triterpenoid core, highly oxygenated, polycyclic and contains a furan ring)"

# Example usage (uncomment the following lines to test):
# test_smiles = "COC(=O)C[C@H]1[C@]2(C)C[C@@]3(O)[C@]1(C)[C@H]1CC[C@@]4(C)[C@@H](OC(=O)[C@H](OC(=O)C(C)(C)O)C4=C1[C@H](OC(C)=O)[C@@]3(O)[C@H]2OC(=O)C(\\C)=C\\C)c1ccoc1"
# result, reason = is_limonoid(test_smiles)
# print(result, reason)