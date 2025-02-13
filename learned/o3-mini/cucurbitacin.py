"""
Classifies: CHEBI:16219 cucurbitacin
"""
"""
Classifies: Cucurbitacin-type compounds
Heuristic improvements:
 - The molecule must parse and should consist only of C, H, and O (cucurbitacins are triterpenoids).
 - It must contain at least 4 rings (a tetracyclic core). Glycosylated variants usually display ≥5 rings.
 - It must have a sufficiently high carbon count (≥25 carbons).
 - Its molecular weight should be above ~300 Da.
 - Typically an α,β–unsaturated ketone (enone; SMARTS "C=CC(=O)") is present. If not, then extra rings (≥5) are required.
 - The molecule must be oxygenated.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_cucurbitacin(smiles: str):
    """
    Determines if a molecule is a cucurbitacin-type compound based on its SMILES string.

    Uses the heuristics:
      - Must be parseable.
      - Only allowed atoms: C (6), H (1), and O (8) (this excludes many false positives).
      - Has at least 4 rings (≥4) (if lacking a free enone, then ≥5 rings are required).
      - Has a minimum carbon count (≥25) and molecular weight (≥300 Da).
      - Many cucurbitacins have a conjugated enone motif (SMARTS: "C=CC(=O)").
      - The molecule must be oxygenated.
      
    Args:
        smiles (str): Input SMILES string.
    
    Returns:
        (bool, str): True with reason if the criteria pass; else False with an explanation.
    """
    # Parse the SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Reject molecules that contain atoms other than C, H, or O.
    allowed_atomic_nums = {1, 6, 8}
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in allowed_atomic_nums:
            return False, (f"Contains heteroatom {atom.GetSymbol()} (atomic number {atom.GetAtomicNum()}) "
                           "which is not typical for cucurbitacin-type compounds")
    
    # Count rings using the RDKit ring info.
    ring_info = mol.GetRingInfo()
    total_rings = ring_info.NumRings()
    if total_rings < 4:
        return False, f"Only {total_rings} rings detected; cucurbitacins typically are tetracyclic (>=4 rings)"
    
    # Count the number of carbon atoms.
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 25:
        return False, f"Only {carbon_count} carbon atoms detected; too few for a cucurbitacin skeleton or derivative"
    
    # Calculate the molecular weight.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:
        return False, f"Molecular weight ({mol_wt:.1f} Da) is too low for a cucurbitacin derivative"
    
    # Check for presence of oxygen atoms.
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if oxygen_count < 1:
        return False, "No oxygen atoms detected; cucurbitacins are oxygenated compounds"
    
    # Look for the conjugated enone motif.
    # This SMARTS looks for a C=C-C(=O) fragment.
    enone_smarts = "C=CC(=O)"
    enone_pattern = Chem.MolFromSmarts(enone_smarts)
    enone_found = mol.HasSubstructMatch(enone_pattern)
    
    # If no free enone motif is found, require extra rings (≥5) to allow for glycosylated variants.
    if not enone_found and total_rings < 5:
        return False, "No conjugated enone motif detected and insufficient rings (<5) to suggest a glycosylated cucurbitacin"
    
    status = "contains a conjugated enone motif" if enone_found else "lacks a free enone but appears glycosylated (>=5 rings)"
    reason = (f"Molecule has {total_rings} rings, {carbon_count} carbons, molecular weight {mol_wt:.1f} Da, "
              f"{status}, and is oxygenated, consistent with cucurbitacin-type compounds.")
    
    return True, reason

# Example usage (for testing):
if __name__ == '__main__':
    # Test with cucurbitacin I (non-glycosylated) SMILES
    test_smiles = "CC(C)(O)\\C=C\\C(=O)[C@](C)(O)[C@H]1[C@H](O)C[C@@]2(C)[C@@H]3CC=C4[C@@H](C=C(O)C(=O)C4(C)C)[C@]3(C)C(=O)C[C@]12C"
    result, explanation = is_cucurbitacin(test_smiles)
    print(result)
    print(explanation)