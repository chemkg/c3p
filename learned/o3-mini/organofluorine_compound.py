"""
Classifies: CHEBI:37143 organofluorine compound
"""
"""
Classifies: Organofluorine Compound
An organofluorine compound is defined as any compound containing at least one carbon–fluorine bond.
This improved version uses a SMARTS pattern "[#6]-[F]" to look for an explicit C–F bond.
"""
from rdkit import Chem

def is_organofluorine_compound(smiles: str):
    """
    Determines if a molecule is an organofluorine compound based on its SMILES string.
    An organofluorine compound is defined as any compound that contains at least one carbon–fluorine (C–F) bond.
    
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
    
    # Define a SMARTS pattern for a carbon-fluorine bond.
    # In SMARTS, [#6] represents any carbon atom and [F] a fluorine atom.
    cf_pattern = Chem.MolFromSmarts("[#6]-[F]")
    if cf_pattern is None:
        return False, "Error in SMARTS pattern"
    
    # Check if the molecule contains at least one match to the pattern.
    if mol.HasSubstructMatch(cf_pattern):
        return True, "Contains at least one carbon–fluorine bond"
    else:
        return False, "No carbon–fluorine bonds detected"

# Example usage (for local testing; uncomment to run):
# test_smiles = [
#     "C[Si](Cn1cncn1)(c1ccc(F)cc1)c1ccc(F)cc1",  # flusilazole; expected True
#     "CCC",  # no C–F bond; expected False
# ]
# for smi in test_smiles:
#     result, reason = is_organofluorine_compound(smi)
#     print(f"SMILES: {smi}\nResult: {result}\nReason: {reason}\n")