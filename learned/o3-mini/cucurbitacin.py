"""
Classifies: CHEBI:16219 cucurbitacin
"""
"""
Classifies: Cucurbitacin-type compounds (tetracyclic triterpenoids derived from cucurbitane)

Heuristic notes (updated):
1. Cucurbitacins have a tetracyclic (≥4 rings) core. Many glycosylated derivatives have extra rings (≥5).
2. They have a high carbon count; typical aglycones have around 30 carbons. Including sugars can raise that number.
3. They are oxygenated.
4. Many (but not all) cucurbitacins display an enone motif (an α,β-unsaturated carbonyl, as in C=CC(=O)).
5. Because glycosylation may “mask” the free enone motif, if the enone is not found we require additional rings (≥5)
   and a higher carbon count (≥25) to capture plausible glycosides while excluding most false positives.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_cucurbitacin(smiles: str):
    """
    Determines if a molecule is a cucurbitacin-type compound based on its SMILES string.
    
    Heuristics used:
      - The molecule must be parsable.
      - It should contain at least 4 rings (a tetracyclic core). Glycosylated compounds often have ≥5 rings.
      - It should have a sufficiently high carbon count (≥25 atoms) to be consistent with the core cucurbitacin skeleton
        (or its glycosylated variants).
      - Its molecular weight should be above ~300 Da.
      - Many cucurbitacins contain a conjugated enone motif (C=CC(=O)). If the enone motif is not found,
        then we require the molecule to have extra rings (≥5) to allow for glycosylated derivatives.
      - The molecule must be oxygenated.
      
    Returns:
        (bool, str): True with a reason if the molecule passes the criteria; otherwise False with an explanation.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Count rings using RDKit's ring info.
    ring_info = mol.GetRingInfo()
    total_rings = ring_info.NumRings()
    if total_rings < 4:
        return False, f"Only {total_rings} rings detected; cucurbitacins typically are tetracyclic (>=4 rings)"
    
    # Count carbon atoms.
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 25:
        return False, f"Only {carbon_count} carbon atoms detected; too few for a typical cucurbitacin core or derivative"
    
    # Calculate molecular weight.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:
        return False, f"Molecular weight ({mol_wt:.1f} Da) is too low for a cucurbitacin derivative"
    
    # Check for enone motif via SMARTS search.
    # This SMARTS looks for a conjugated α,β-unsaturated carbonyl: C=C-C(=O)
    enone_smarts = "C=CC(=O)"
    enone_pattern = Chem.MolFromSmarts(enone_smarts)
    enone_found = mol.HasSubstructMatch(enone_pattern)
    
    # If no enone motif is found, require extra rings to allow for glycosylated variants
    if not enone_found:
        if total_rings < 5:
            return False, "No conjugated enone motif detected and insufficient rings (<5) to suggest a glycosylated cucurbitacin"
    
    # Check for the presence of oxygen atoms.
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if oxygen_count < 1:
        return False, "No oxygen atoms detected; cucurbitacins are oxygenated compounds"
    
    reason = ("Molecule has {0} rings, {1} carbons, molecular weight {2:.1f} Da, "
              "{3} an enone motif, and is oxygenated, consistent with cucurbitacin-type compounds."
             ).format(
                 total_rings,
                 carbon_count,
                 mol_wt,
                 "contains" if enone_found else "lacks a free enone but appears to be glycosylated (≥5 rings)"
             )
    
    return True, reason

# Example usage (for testing):
if __name__ == '__main__':
    # Test with cucurbitacin I (non-glycosylated) SMILES
    test_smiles = "CC(C)(O)\\C=C\\C(=O)[C@](C)(O)[C@H]1[C@H](O)C[C@@]2(C)[C@@H]3CC=C4[C@@H](C=C(O)C(=O)C4(C)C)[C@]3(C)C(=O)C[C@]12C"
    result, explanation = is_cucurbitacin(test_smiles)
    print(result)
    print(explanation)