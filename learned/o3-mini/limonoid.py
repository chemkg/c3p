"""
Classifies: CHEBI:39434 limonoid
"""
"""
Classifies: Limonoid
Definition:
  Any triterpenoid that is highly oxygenated and has a prototypical structure either containing or derived from a precursor 
  with a 4,4,8-trimethyl-17-furanylsteroid skeleton. (e.g. limonin, azadirachtin, and other structurally related compounds).
  
This implementation uses heuristic criteria:
  - The molecule must be in a reasonable molecular weight range for a tetranortriterpenoid (300-700 Da)
  - Carbon count is expected to be between ~25 and 40 (typical for a triterpenoid skeleton).
  - The molecule must be oxygen-rich (at least 4 oxygen atoms), and if oxygenation is moderate (4-6), it must contain
    a furan ring (SMARTS "c1ccoc1"). If oxygenation is high (>=7) the explicit furan match is not strictly required.
  - The molecule must be polycyclic (at least 4 rings).
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_limonoid(smiles: str):
    """
    Determines if a molecule is a limonoid based on its SMILES string.
    
    Heuristic strategy:
      - Parse the SMILES.
      - Check that the molecular weight is between 300 and 700 Da.
      - Verify that the number of carbon atoms is roughly between 25 and 40.
      - Check that the molecule is highly oxygenated by requiring at least 4 oxygen atoms.
      - If oxygen count is moderate (4 to 6), require a furan ring substructure (SMARTS "c1ccoc1").
        Otherwise, if oxygenation is very high (>=7), the absence of a furan ring may be allowed.
      - Ensure the molecule has a polycyclic core (at least 4 rings).
      
    Args:
       smiles (str): SMILES representation of the molecule.
       
    Returns:
       bool: True if the molecule is classified as a limonoid, otherwise False.
       str: Explanation of the classification decision.
    """
    
    # Parse SMILES into an rdkit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Calculate molecular weight and require it to be within a range typical for limonoids.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300 or mol_wt > 700:
        return False, f"Molecular weight ({mol_wt:.1f} Da) is outside the expected range (300-700 Da)"
    
    # Count the number of carbon atoms (limonoids are derived from triterpenoids, usually 25-40 carbons)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 25 or c_count > 40:
        return False, f"Carbon count ({c_count}) is outside the typical range for a triterpenoid skeleton (25-40)"
    
    # Count oxygen atoms; limonoids tend to be highly oxygenated.
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 4:
        return False, f"Not enough oxygen atoms ({o_count}); expected at least 4 for a limonoid"
    
    # Check for polycyclic structure: require at least 4 rings.
    n_rings = mol.GetRingInfo().NumRings()
    if n_rings < 4:
        return False, f"Not enough rings ({n_rings}); expected a polycyclic structure typical of limonoids"
    
    # Check for the presence of a furan ring using a SMARTS pattern.
    furan_pattern = Chem.MolFromSmarts("c1ccoc1")
    has_furan = mol.HasSubstructMatch(furan_pattern)
    
    # If oxygen count is moderate, then an explicit furan ring is required.
    if 4 <= o_count <= 6 and not has_furan:
        return False, "Furan ring not found (and oxygenation is only moderate); expected typical furan ring in limonoids"
    
    # Otherwise, if oxygenation is high, we allow the absence of a clear furan ring.
    # (This accounts for limonoids derived from a furanylsteroid precursor that have lost the explicit furan moiety.)
    
    return True, "Molecule appears to be a limonoid (triterpenoid core, highly oxygenated, polycyclic and (if oxygenation is moderate) contains a furan ring)"

# Example usage (uncomment to test):
# test_smiles = "COC(=O)C[C@H]1[C@]2(C)C[C@@]3(O)[C@]1(C)[C@H]1CC[C@@]4(C)[C@@H](OC(=O)[C@H](OC(=O)C(C)(C)O)C4=C1[C@H](OC(C)=O)[C@@]3(O)[C@H]2OC(=O)C(\\C)=C\\C)c1ccoc1"
# result, reason = is_limonoid(test_smiles)
# print(result, reason)