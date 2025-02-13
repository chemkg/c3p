"""
Classifies: CHEBI:26199 polyprenol
"""
#!/usr/bin/env python
"""
Classifies: polyprenol – any member of the class of prenols possessing the general formula
H-[CH2C(Me)=CHCH2]nOH (with n >= 2, i.e. more than one isoprene unit).
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_polyprenol(smiles: str):
    """
    Determines if a molecule is a polyprenol based on its SMILES string.
    Polyprenols are (typically linear) isoprenoids that consist of more than one isoprene unit [CH2C(Me)=CHCH2]
    and contain at least one hydroxyl group.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is classified as a polyprenol, False otherwise.
        str: A textual reason for the classification decision.
    """
    # Parse SMILES into a molecule; if parsing fails, return error.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check that the molecule has at least one hydroxyl group.
    # The SMARTS "[OX2H]" matches an oxygen with two connections and a hydrogen, as found in -OH.
    hydroxyl_pat = Chem.MolFromSmarts("[OX2H]")
    if not mol.HasSubstructMatch(hydroxyl_pat):
        return False, "No hydroxyl (-OH) group found"

    # It is expected that polyprenol structures are built from multiple isoprene units.
    # The ideal repeating unit is: CH2-C(Me)=CH-CH2.
    # However, at the chain ends the pattern is modified.
    # We define three SMARTS patterns:
    #   1. pattern_internal: the complete unit that can appear in the middle: [CH2]C([CH3])=C[CH2]
    #   2. pattern_terminal_omega: a terminal unit at the omega end: [CH2]C([CH3])=C[CH3]
    #   3. pattern_terminal_alpha: a terminal unit at the alpha end (adjacent to -OH): [OH]CC([CH3])=C
    pattern_internal = Chem.MolFromSmarts("[CH2]C([CH3])=C[CH2]")
    pattern_term_omega = Chem.MolFromSmarts("[CH2]C([CH3])=C[CH3]")
    pattern_term_alpha = Chem.MolFromSmarts("[OH]CC([CH3])=C")

    matches_internal = mol.GetSubstructMatches(pattern_internal)
    matches_term_omega = mol.GetSubstructMatches(pattern_term_omega)
    matches_term_alpha = mol.GetSubstructMatches(pattern_term_alpha)
    
    # Count the total number of detected isoprene‐like units.
    total_units = len(matches_internal) + len(matches_term_omega) + len(matches_term_alpha)
    
    # To be considered a polyprenol the molecule should have more than one isoprene unit.
    if total_units < 2:
        return False, f"Only {total_units} isoprene unit(s) detected (need at least 2)"
    
    # Optional: Check that the molecule has a sufficiently long carbon chain.
    # The general formula H-[CH2C(Me)=CHCH2]nOH contains 5*n+2 carbons.
    # We can require that number of carbon atoms is at least 10.
    c_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    if len(c_atoms) < 10:
        return False, f"Too few carbons ({len(c_atoms)}) for a polyprenol structure"

    # You may also check other features (e.g., molecular weight) if desired.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 150:
        return False, f"Molecular weight ({mol_wt:.1f} Da) is too low for a polyprenol"

    return True, f"Polyprenol detected with {total_units} isoprene-like unit(s) and at least one -OH group."

# Example usage:
if __name__ == '__main__':
    # Test with a few SMILES examples.
    test_smiles = {
        "(E,E,E)-geranylgeraniol": "CC(C)=CCC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C\\CO",
        "geraniol": "CC(C)=CCC\\C(C)=C\\CO",
        "(2E,6E,10E)-omega-hydroxyfarnesol": "C\\C(CO)=C/CC\\C(C)=C\\CC\\C(C)=C\\CO"
    }
    for name, s in test_smiles.items():
        result, reason = is_polyprenol(s)
        print(f"{name}: {result} ({reason})")