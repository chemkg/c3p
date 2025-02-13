"""
Classifies: CHEBI:26244 prenols
"""
#!/usr/bin/env python3
"""
Classifies: Prenols
Definition: Any alcohol possessing the general formula H-[CH2C(Me)=CHCH2]nOH in which the carbon
skeleton is composed of one or more isoprene units. Prenol derivatives (such as phosphorylated species)
are accepted provided the terminal alcohol (free or phosphorylated) is present and the overall
carbon count matches that expected for a prenol (i.e. 5 * n carbons, with n>=1).
"""

from rdkit import Chem

def is_prenols(smiles: str):
    """
    Determines if a molecule belongs to the prenol class based on its SMILES string.
    
    Classification strategy:
      1. Parse the SMILES string.
      2. Count carbons and require that the number of carbons is a multiple of 5 (except for prenol, 5 C, which is allowed).
      3. Look for a terminal alcohol function. We require that the terminal –OH is attached to a CH2.
         We also accept a phosphorylated derivative (e.g. CH2O–P...) by a separate SMARTS.
      4. For molecules larger than prenol (i.e. with 10 or more carbons) we search for at least one full
         isoprene repeating unit (CH2–C(CH3)=CH–CH2) using an appropriate SMARTS. Since our expected
         count of repeating units (as indicated in our example outcomes) is defined as:
              expected_units = 1 for prenol (5 carbons)
              expected_units = (total_C/5) - 1   for molecules with total_C >=10.
         We require that the number of non-overlapping isoprene matches is at least expected_units.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule is classified as a prenol (or prenol derivative), False otherwise.
        str: A brief reason describing the classification result.
    """
    # Convert SMILES to molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count carbon atoms.
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)

    # Check for terminal alcohol. We require that the hydroxyl is on a CH2.
    # Pattern 1: free terminal alcohol: [CH2][OX2H]
    terminal_alc = Chem.MolFromSmarts("[CH2][OX2H]")
    # Pattern 2: phosphorylated terminal alcohol: [CH2]O[P]
    terminal_alc_phos = Chem.MolFromSmarts("[CH2]O[P]")
    
    has_terminal_alcohol = mol.HasSubstructMatch(terminal_alc) or mol.HasSubstructMatch(terminal_alc_phos)
    if not has_terminal_alcohol:
        return False, "No terminal alcohol (free or phosphorylated) found"

    # For prenol (the smallest case, 5 carbons) we accept the molecule, even though the isoprene unit SMARTS
    # may not fire.
    if c_count == 5:
        return True, "Found terminal alcohol and 1 isoprene repeating unit(s)."
    
    # For molecules with 10 or more carbons, enforce that the overall carbon count is a multiple of 5.
    if c_count % 5 != 0:
        return False, f"Carbon count {c_count} is not a multiple of 5, which is atypical for prenols."

    # Determine expected number of repeating isoprene units.
    # Based on the provided examples, a molecule with 10 carbons (geraniol) should have 1 unit,
    # with 15 carbons (farnesol) 2 units, with 20 carbons 3 units, etc.
    expected_units = (c_count // 5) - 1

    # For molecules with 10 or more carbons, check for the isoprene repeating unit.
    # We use a SMARTS pattern that represents a full isoprene unit: CH2–C(CH3)=CH–CH2.
    # This pattern is defined as: "[CH2]-[C](C)=C-[CH2]". 
    isoprene_pattern = Chem.MolFromSmarts("[CH2]-[C](C)=C-[CH2]")
    matches = mol.GetSubstructMatches(isoprene_pattern)
    
    # Note: in larger systems these matches may overlap, but we require at least as many matches
    # as expected repeating units.
    if len(matches) < expected_units:
        return False, (f"Found terminal alcohol but only {len(matches)} isoprene repeating unit(s) detected; "
                       f"expected at least {expected_units} based on carbon count {c_count}.")
    
    return True, f"Found terminal alcohol and {expected_units} isoprene repeating unit(s)."

# When run as a script, execute a few test cases.
if __name__ == '__main__':
    test_smiles_list = [
        "CC(C)=CCC\\C(C)=C\\CO",  # (E,E,E)-geranylgeraniol (20 C, expect 3 units)
        "C\\C(CO)=C/CC\\C(C)=C\\CC\\C(C)=C\\CO",  # (2E,6E,10E)-omega-hydroxyfarnesol (15 C, expect 2 units)
        "CC(C)=CCC\\C(C)=C\\CO",  # geraniol (10 C, expect 1 unit)
        "CC(C)=CCO",  # prenol (5 C, expect 1 unit)
    ]
    for smi in test_smiles_list:
        result, reason = is_prenols(smi)
        print(f"SMILES: {smi}\n  Is prenol? {result}\n  Reason: {reason}\n")