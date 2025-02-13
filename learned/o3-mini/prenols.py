"""
Classifies: CHEBI:26244 prenols
"""
#!/usr/bin/env python3
"""
Classifies: Prenols
Definition: Any alcohol possessing the general formula H-[CH2C(Me)=CHCH2]nOH in which the carbon
skeleton is composed of one or more isoprene units. Prenol derivatives (e.g. phosphorylated species)
are accepted provided the terminal alcohol (free or phosphorylated) is present and the overall
carbon count is consistent with a prenol derivative (i.e. 5*n carbons, with n>=1).

Revision notes:
  • In the previous approach, the number of isoprene substructure matches (using a SMARTS pattern)
    was used to “count” repeating units. However, for larger and highly branched prenols these matches
    is often under‐detected.
  • To improve the classification we now require:
      1. A valid molecule that is acyclic (prenols are linear polyprenyl chains).
      2. A total carbon count that is either 5 or (if >5) a multiple of 5.
      3. A terminal alcohol is present – either free ([CH2][OX2H]) or phosphorylated ([CH2]O[P]).
      4. (For molecules larger than prenol) the presence of at least one carbon–carbon double bond
         and at least one methyl‐substituted double bond ([C]([CH3])=[C]) is required.
  • The expected number of isoprene units is then computed as:
          expected_units = 1    if total carbons == 5
          expected_units = (c_count // 5) - 1  if c_count >= 10
    and this is reported.
"""

from rdkit import Chem

def is_prenols(smiles: str):
    """
    Determines if a molecule belongs to the prenol class (or prenol derivative) based on its SMILES.
    
    The algorithm:
      1. Parse the SMILES string.
      2. Reject if the molecule contains rings (prenol chains are acyclic).
      3. Count carbon atoms. For prenols, the total carbon count should be exactly 5 (for prenol)
         or a multiple of 5 (for longer prenols).
      4. Check for the presence of a terminal alcohol group, either free or phosphorylated.
         A terminal –OH is defined as a hydroxyl attached to a CH2 group.
      5. For molecules larger than prenol (i.e. >5 carbons), also require that at least one C=C
         bond exists and that a methyl-substituted double bond pattern ([C]([CH3])=[C]) is present.
         This is a diagnostic feature of the isoprene unit.
      6. Compute the expected number of isoprene repeating units:
             if total_carbons == 5: (prenol) report 1 unit,
             else expected_units = (total_carbons // 5) - 1.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if classified as a prenol (or derivative), False otherwise.
        str: A reason message describing the classification result.
    """

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Reject molecules with rings as prenols are linear chains.
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule contains rings which is atypical for prenols"
    
    # Count carbon atoms.
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 5:
        return False, f"Too few carbons ({c_count}) to be a prenol"
    
    # For molecules > prenol, the carbon count must be a multiple of 5.
    if c_count > 5 and (c_count % 5 != 0):
        return False, f"Carbon count {c_count} is not a multiple of 5, which is atypical for prenols"
    
    # Check for terminal alcohol: free -OH ([CH2][OX2H]) or phosphorylated ([CH2]O[P])
    terminal_alc = Chem.MolFromSmarts("[CH2][OX2H]")
    terminal_alc_phos = Chem.MolFromSmarts("[CH2]O[P]")
    has_terminal_alcohol = mol.HasSubstructMatch(terminal_alc) or mol.HasSubstructMatch(terminal_alc_phos)
    if not has_terminal_alcohol:
        return False, "No terminal alcohol (free or phosphorylated) found"
    
    # For molecules with >5 carbons (i.e. longer prenols), check for isoprenoid features:
    # Require presence of a carbon–carbon double bond.
    dbl_bond = Chem.MolFromSmarts("[C]=[C]")
    if c_count > 5 and not mol.HasSubstructMatch(dbl_bond):
        return False, "No C=C double bond found, required for an isoprenoid structure"
        
    # Also require a methyl-substituted double bond, e.g. [C]([CH3])=[C]
    methyl_dbl = Chem.MolFromSmarts("[C]([CH3])=[C]")
    if c_count > 5 and not mol.HasSubstructMatch(methyl_dbl):
        return False, "No methyl-substituted double bond ([C]([CH3])=[C]) detected, not typical of isoprene units"
    
    # Determine the expected number of isoprene units.
    if c_count == 5:
        expected_units = 1
    else:
        expected_units = (c_count // 5) - 1

    return True, f"Found terminal alcohol; acyclic molecule with {c_count} carbons, consistent with {expected_units} isoprene repeating unit(s)."

# When run as a script, we may test a few examples.
if __name__ == '__main__':
    test_smiles_list = [
        "CC(C)=CCC\\C(C)=C\\CO",   # (E,E,E)-geranylgeraniol, 20 C => expected 3 repeating units.
        "C\\C(CO)=C/CC\\C(C)=C\\CC\\C(C)=C\\CO",   # (2E,6E,10E)-omega-hydroxyfarnesol, 15 C => expected 2 repeating units.
        "CC(C)=CCC\\C(C)=C\\CO",  # geraniol, 10 C => expected 1 repeating unit.
        "CC(C)=CCO",  # prenol, 5 C => expected 1 unit.
        "CC(=O)C(O)C(=O)COP([O-])([O-])=O",  # false positive from previous approach.
    ]
    for smi in test_smiles_list:
        result, reason = is_prenols(smi)
        print(f"SMILES: {smi}\n  Is prenol? {result}\n  Reason: {reason}\n")