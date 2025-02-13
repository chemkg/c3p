"""
Classifies: CHEBI:26244 prenols
"""
#!/usr/bin/env python3
"""
Classifies: Prenols
Definition: Any alcohol possessing the general formula H-[CH2C(Me)=CHCH2]nOH in which the carbon
skeleton is composed of one or more isoprene units (biogenetic precursors of the isoprenoids).
Some prenol derivatives exist as phosphorylated species (e.g. diphosphates). This function
checks for both a terminal alcohol (free or phosphorylated) and the presence of one or more
isoprene repeating units.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_prenols(smiles: str):
    """
    Determines if a molecule belongs to the prenols class based on its SMILES string.
    The approach is:
      1. Parse the SMILES string.
      2. Verify that the molecule has a terminal alcohol feature. Either a free hydroxyl group
         (–OH) attached to a carbon, or an oxygen attached to a phosphate group (indicating a
         phosphorylated alcohol derivative) is accepted.
      3. Look for at least one occurrence of the isoprene repeating unit:
             CH2–C(CH3)=CH–CH2
         using the SMARTS "[CH2]-[C]([CH3])=C-[CH2]".
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule is classified as a prenol (or a prenol derivative), False otherwise.
        str: A brief reason describing the classification result.
    """
    # Parse SMILES into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of an alcohol group.
    # We define two patterns: one for a free terminal alcohol group (e.g. CH2OH)
    # and one for a phosphorylated derivative (e.g. CH2O–P(=O)(O)...) which arises from phosphorylation.
    alcohol_pattern = Chem.MolFromSmarts("[#6][OX2H]")  # A carbon atom directly bound to a hydroxyl group
    phospho_pattern = Chem.MolFromSmarts("[#6]O[P]")     # A carbon connected to an oxygen then to phosphorus
    
    has_alcohol = mol.HasSubstructMatch(alcohol_pattern) or mol.HasSubstructMatch(phospho_pattern)
    if not has_alcohol:
        return False, "No terminal alcohol (or phosphorylated alcohol) moiety found"

    # Check for one or more isoprene repeating units.
    # The pattern "[CH2]-[C]([CH3])=C-[CH2]" represents CH2 – C(CH3)=CH – CH2
    isoprene_pattern = Chem.MolFromSmarts("[CH2]-[C]([CH3])=C-[CH2]")
    isoprene_matches = mol.GetSubstructMatches(isoprene_pattern)
    if len(isoprene_matches) < 1:
        return False, "No isoprene repeating unit (CH2–C(CH3)=CH–CH2) found"

    # As an extra check, we might want to verify that the number of carbons is in line with a chain
    # composed of 5-carbon units. The general formula requires 5*n carbon atoms. We warn if this is unusual.
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count % 5 != 0:
        # This is a warning rather than a rejection since natural products may have additional modifications.
        return True, f"Found terminal alcohol and isoprene unit(s), but carbon count {c_count} is not a typical multiple of 5."

    return True, f"Found terminal alcohol and {len(isoprene_matches)} isoprene repeating unit(s)."

# If the file is run as a script, one might include some test calls.
if __name__ == "__main__":
    # Test examples (you can add more to this list)
    test_smiles_list = [
        "CC(C)=CCC\\C(C)=C\\CO",  # geranylgeraniol analogue
        "CC(C)=CCO",            # prenol example
        "C1CCCCC1",             # cyclohexane (should fail)
    ]
    for smi in test_smiles_list:
        result, reason = is_prenols(smi)
        print(f"SMILES: {smi}\n  Is prenol? {result}\n  Reason: {reason}\n")