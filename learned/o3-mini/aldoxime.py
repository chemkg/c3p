"""
Classifies: CHEBI:22307 aldoxime
"""
"""
Classifies: Aldoxime (Oximes of aldehydes RCH=NOH)

An aldoxime is an oxime derivative of an aldehyde,
meaning the oxime group is attached to a carbon that once was an aldehyde.
For an aldehyde, the carbon has exactly one hydrogen (RCH=O goes to RCH=NOH).
This script first searches for the general C=N–OH pattern,
then verifies that the carbon atom involved has exactly one hydrogen.
"""

from rdkit import Chem

def is_aldoxime(smiles: str):
    """
    Determines if a molecule is an aldoxime based on its SMILES string.
    An aldoxime is characterized by an oxime group derived from an aldehyde,
    i.e., it must have the functional group RCH=NOH where the carbon atom has exactly one hydrogen.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is classified as an aldoxime, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string into an RDKit molecule object.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens so that hydrogen counts are correctly available.
    mol = Chem.AddHs(mol)
    
    # Define a general SMARTS pattern for an oxime group: a carbon double-bonded to a nitrogen
    # that is bonded to a hydroxyl group.
    # We do not include explicit hydrogen count on the carbon in the pattern.
    pattern = Chem.MolFromSmarts("[#6]=[NX2][OX1H]")
    if pattern is None:
        return False, "Error in SMARTS pattern"
    
    matches = mol.GetSubstructMatches(pattern)
    if not matches:
        return False, "Aldoxime group not found in the molecule"
    
    # For each matched pattern, verify that the carbon (first atom in match) has exactly one hydrogen.
    for match in matches:
        carbon_idx = match[0]
        carbon_atom = mol.GetAtomWithIdx(carbon_idx)
        # GetTotalNumHs returns the total number of hydrogens (implicit + explicit).
        if carbon_atom.GetTotalNumHs() == 1:
            return True, "Contains an aldoxime group (RCH=NOH)"
    
    # If no match has a carbon with exactly one hydrogen attached, return false.
    return False, "Aldoxime group found but the aldehyde carbon does not have exactly one hydrogen"

# Uncomment below lines for simple testing examples:
# test_smiles = ["[H]\\C(C(C)C)=N/O", "C([C@@H](/C(=N/O)/[H])C)C", "COc1cc(\\C=N/O)nc(c1)-c1ccccn1"]
# for smi in test_smiles:
#     result, reason = is_aldoxime(smi)
#     print(f"SMILES: {smi}\nResult: {result}\nReason: {reason}\n")