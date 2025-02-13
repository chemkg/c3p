"""
Classifies: CHEBI:37143 organofluorine compound
"""
"""
Classifies: Organofluorine Compound
An organofluorine compound is defined as any compound that contains at least one carbon–fluorine (C–F) bond.
This version directly inspects all bonds in the molecule.
"""
from rdkit import Chem

def is_organofluorine_compound(smiles: str):
    """
    Determines if a molecule is an organofluorine compound based on its SMILES string.
    An organofluorine compound is defined as any organic compound that contains at least one carbon–fluorine (C–F) bond.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule contains at least one C–F bond, False otherwise.
        str: Reason for classification.
    """
    # Parse the SMILES string into an RDKit molecule object.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Iterate over every bond in the molecule.
    # For each bond, check if one atom is carbon (atomic number 6)
    # and the other is fluorine (atomic number 9). This method does not depend
    # on the bond type denoted in the SMARTS pattern and is more robust.
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        if (a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 9) or (a1.GetAtomicNum() == 9 and a2.GetAtomicNum() == 6):
            return True, "Contains at least one carbon–fluorine bond"
    
    return False, "No carbon–fluorine bonds detected"

# Example usage (for local testing):
# test_smiles = [
#     "C[Si](Cn1cncn1)(c1ccc(F)cc1)c1ccc(F)cc1",  # flusilazole; expected True
#     "CCC",  # no C–F bond; expected False
# ]
# for smi in test_smiles:
#     result, reason = is_organofluorine_compound(smi)
#     print(f"SMILES: {smi}\nResult: {result}\nReason: {reason}\n")