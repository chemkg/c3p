"""
Classifies: CHEBI:37143 organofluorine compound
"""
"""
Classifies: Organofluorine Compound
An organofluorine compound is defined as any compound containing at least one carbon–fluorine bond.
In this improved version we inspect every bond to ensure that a C–F bond is found.
"""
from rdkit import Chem

def is_organofluorine_compound(smiles: str):
    """
    Determines if a molecule is an organofluorine compound based on its SMILES string.
    An organofluorine compound is defined (for our purposes) as any compound that contains
    at least one carbon–fluorine (C–F) bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule contains at least one C–F bond, False otherwise.
        str: Reason for classification.
    """
    # Parse the SMILES string into a molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Iterate over each bond and check if one atom is carbon (atomic num 6) and the other is fluorine (atomic num 9)
    for bond in mol.GetBonds():
        atom1 = bond.GetBeginAtom()
        atom2 = bond.GetEndAtom()
        # Check both directions
        if (atom1.GetAtomicNum() == 6 and atom2.GetAtomicNum() == 9) or (atom1.GetAtomicNum() == 9 and atom2.GetAtomicNum() == 6):
            return True, "Contains at least one carbon–fluorine bond"
    
    # If we found no such bond, then the compound is not organofluorine.
    return False, "No carbon–fluorine bonds detected"

# Example usage (uncomment for testing):
# test_smiles = [
#     "C[Si](Cn1cncn1)(c1ccc(F)cc1)c1ccc(F)cc1",  # flusilazole; expected True
#     "CCC",  # no C–F bond; expected False
# ]
# for smi in test_smiles:
#     result, reason = is_organofluorine_compound(smi)
#     print(f"SMILES: {smi}\nResult: {result}\nReason: {reason}\n")