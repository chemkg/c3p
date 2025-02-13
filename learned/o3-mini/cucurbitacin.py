"""
Classifies: CHEBI:16219 cucurbitacin
"""
"""
Classifies: Cucurbitacin-type compounds (tetracyclic triterpenoids derived from cucurbitane)

Heuristic notes:
1. Cucurbitacins are tetracyclic (contain at least four rings, though glycosylated derivatives can add extra rings).
2. They typically have a fairly high carbon count (usually around 30 on the core skeleton, though extra groups may add more).
3. They tend to be oxygenated (one or more oxygen atoms) and display a conjugated enone motif (i.e. a C=C–C(=O) fragment).
4. These criteria are heuristic and may not cover all cucurbitacin derivatives.

The function returns (True, reason) if all criteria are met and (False, reason) otherwise.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_cucurbitacin(smiles: str):
    """
    Determines if a molecule is a cucurbitacin-type compound based on its SMILES string.
    
    Cucurbitacins are tetracyclic triterpenoids derived from cucurbitane.
    Heuristics used:
      - Molecule must be parsed successfully.
      - The molecule should contain at least 4 rings.
      - The carbon count should be high (heuristically >=20 to allow for the 30-carbon core or glycosylated derivatives).
      - Molecular weight should be above ~300 Da.
      - An enone motif (C=C–C(=O)) must be present, which is a frequently observed feature of these molecules.
      - At least one oxygen atom is required.
      
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as cucurbitacin, False otherwise.
        str: Explanation of the classification decision.
    """
    # Parse the SMILES string into a molecule object.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Get ring information.
    ring_info = mol.GetRingInfo()
    total_rings = ring_info.NumRings()
    if total_rings < 4:
        return False, f"Only {total_rings} rings detected; cucurbitacins are tetracyclic (>=4 rings)"
    
    # Count carbon atoms in the molecule.
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 20:
        return False, f"Only {carbon_count} carbon atoms detected; too few for a cucurbitacin core"
    
    # Calculate molecular weight.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:
        return False, f"Molecular weight ({mol_wt:.1f} Da) is too low for cucurbitacin derivatives"
    
    # Look for an enone motif (a C=C-C(=O) fragment), frequently found in cucurbitacins.
    enone_pattern = Chem.MolFromSmarts("C=C-C(=O)")
    if not mol.HasSubstructMatch(enone_pattern):
        return False, "No enone motif (C=C-C(=O)) detected, which is common in cucurbitacins"
    
    # Count oxygen atoms.
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if oxygen_count < 1:
        return False, "No oxygen atoms detected; cucurbitacins are oxygenated compounds"
    
    # If all heuristics are passed, we classify this molecule as a cucurbitacin derivative.
    reason = ("Molecule has {0} rings, {1} carbons, molecular weight {2:.1f} Da, contains an enone motif, "
              "and is oxygenated, consistent with cucurbitacin-type compounds."
             ).format(total_rings, carbon_count, mol_wt)
    
    return True, reason

# Example usage:
if __name__ == '__main__':
    # SMILES for cucurbitacin I as an example:
    test_smiles = "CC(C)(O)\\C=C\\C(=O)[C@](C)(O)[C@H]1[C@H](O)C[C@@]2(C)[C@@H]3CC=C4[C@@H](C=C(O)C(=O)C4(C)C)[C@]3(C)C(=O)C[C@]12C"
    result, explanation = is_cucurbitacin(test_smiles)
    print(result)
    print(explanation)