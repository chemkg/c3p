"""
Classifies: CHEBI:37143 organofluorine compound
"""
"""
Classifies: Organofluorine Compound
An organofluorine compound is defined as any compound containing at least one carbon-fluorine bond.
"""
from rdkit import Chem

def is_organofluorine_compound(smiles: str):
    """
    Determines if a molecule is an organofluorine compound based on its SMILES string.
    A compound is considered an organofluorine compound if it contains at least one carbon-fluorine (C–F) bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule contains at least one C-F bond, False otherwise.
        str: Reason for classification.
    """
    # Parse SMILES string into a molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern for a C-F bond (C is atomic number 6, F is atomic number 9)
    cf_pattern = Chem.MolFromSmarts("[#6]-[#9]")
    if cf_pattern is None:
        return False, "Failed to create substructure pattern for C-F bond"
    
    # Check if there is at least one C-F bond in the molecule
    if mol.HasSubstructMatch(cf_pattern):
        return True, "Contains at least one carbon-fluorine bond"
    else:
        return False, "No carbon-fluorine bonds detected"
        
# Example usage (uncomment for testing):
# smiles_list = [
#     "C[Si](Cn1cncn1)(c1ccc(F)cc1)c1ccc(F)cc1",  # flusilazole, an organofluorine compound
#     "CCC",  # example without a C-F bond
# ]
# for smi in smiles_list:
#     result, reason = is_organofluorine_compound(smi)
#     print(f"SMILES: {smi}\nResult: {result}\nReason: {reason}\n")