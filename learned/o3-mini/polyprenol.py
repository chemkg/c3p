"""
Classifies: CHEBI:26199 polyprenol
"""
#!/usr/bin/env python
"""
Classifies: polyprenol – any member of the class of prenols possessing the general formula
H-[CH2C(Me)=CHCH2]nOH where n>=2 (i.e. more than one isoprene unit).
The classifier requires that the molecule is acyclic, contains at least one hydroxyl group,
has a long enough carbon chain, and displays at least two isoprene‐like fragments.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_polyprenol(smiles: str):
    """
    Determines if a molecule is a polyprenol based on its SMILES string.
    A polyprenol should:
      • Contain at least one hydroxyl group (-OH)
      • Be acyclic (linear isoprenoid chain)
      • Have a sufficiently long carbon chain (e.g. at least 10 carbon atoms)
      • Possess at least two repeating isoprene–like units.
      
    We attempt to detect three types of isoprene–like fragments:
      1. An internal unit: [CH2]C([CH3])=C[CH2]
      2. A terminal omega unit: C([CH3])=C[CH2]O   (where –OH is attached at the end)
      3. A terminal alpha unit: [OH]CC([CH3])=C (where –OH is at the beginning)
      
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if classified as a polyprenol, False otherwise.
        str: Textual explanation for the decision.
    """
    # Parse the molecule from SMILES.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # For polyprenols, we expect a linear (acyclic) structure.
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule contains rings; expected a linear isoprenoid chain"

    # Check that the molecule has at least one hydroxyl group.
    hydroxyl_smarts = "[OX2H]"  # matches an -OH group
    hydroxyl_pat = Chem.MolFromSmarts(hydroxyl_smarts)
    if not mol.HasSubstructMatch(hydroxyl_pat):
        return False, "No hydroxyl (-OH) group found"

    # Check that there are a sufficient number of carbon atoms.
    carbons = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    if len(carbons) < 10:
        return False, f"Too few carbons ({len(carbons)}) for a polyprenol structure"
    
    # (Optional) Check molecular weight; polyprenols usually have MW >150 Da.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 150:
        return False, f"Molecular weight ({mol_wt:.1f} Da) is too low for a polyprenol"
    
    # Define SMARTS patterns for the isoprene-like repeating unit.
    # Pattern for an internal unit: ideal repeating unit [CH2]-C(CH3)=CH-[CH2]
    pattern_internal = Chem.MolFromSmarts("[CH2]C([CH3])=C[CH2]")
    # Pattern for a terminal (omega) unit where the hydroxyl appears at the end: C(CH3)=CH-[CH2]O
    pattern_term_omega = Chem.MolFromSmarts("C([CH3])=C[CH2]O")
    # Pattern for a terminal (alpha) unit where the hydroxyl is at the beginning: [OH]CC(CH3)=C
    pattern_term_alpha = Chem.MolFromSmarts("[OH]CC([CH3])=C")
    
    # Count matches from each pattern.
    matches_internal = mol.GetSubstructMatches(pattern_internal)
    matches_term_omega = mol.GetSubstructMatches(pattern_term_omega)
    matches_term_alpha = mol.GetSubstructMatches(pattern_term_alpha)
    
    # Sum total number of detected isoprene-like units.
    total_units = len(matches_internal) + len(matches_term_omega) + len(matches_term_alpha)
    
    # Debug (optional): uncomment to print pattern hits and their counts.
    # print("Internal hits:", matches_internal)
    # print("Terminal omega hits:", matches_term_omega)
    # print("Terminal alpha hits:", matches_term_alpha)
    
    # To be considered a polyprenol the molecule should have more than one isoprene unit.
    if total_units < 2:
        return False, f"Only {total_units} isoprene unit(s) detected (need at least 2)"
    
    # If all tests pass, we assume it is a polyprenol.
    return True, f"Polyprenol detected with {total_units} isoprene-like unit(s) and at least one -OH group."

# Example usage:
if __name__ == "__main__":
    # Test with several example SMILES:
    test_smiles = {
        "(E,E,E)-geranylgeraniol": "CC(C)=CCC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C\\CO",
        "(2E,6E,10E)-omega-hydroxyfarnesol": "C\\C(CO)=C/CC\\C(C)=C\\CC\\C(C)=C\\CO",
        "Bionectin F": "OC(CCC=C(C)C)(CCCC(O)(CCCC(O)(CCCC(O)(CCCC(O)(CCCC(O)(CCCC(O)(CCC/C(=C/CC/C(=C/CO)/C)/C)C)C)C)C)C)C)C",
        "geraniol": "CC(C)=CCC\\C(C)=C\\CO",
        "nerol": "C(=C\\CO)(\\CCC=C(C)C)/C"
    }
    for name, s in test_smiles.items():
        result, reason = is_polyprenol(s)
        print(f"{name}: {result} ({reason})")