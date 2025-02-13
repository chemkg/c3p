"""
Classifies: CHEBI:17792 organohalogen compound
"""
"""
Classifies: Organohalogen compound
Definition: A compound containing at least one carbon-halogen bond (halogen = F, Cl, Br, I).
"""
from rdkit import Chem

def is_organohalogen_compound(smiles: str):
    """
    Determines if a molecule is an organohalogen compound based on its SMILES string.
    The classification is made if the molecule contains at least one carbon-halogen bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule contains at least one carbon-halogen bond, False otherwise.
        str: Reason for the classification.
    """
    # Parse SMILES string into a molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a SMARTS pattern that matches a bond between carbon and any halogen (F, Cl, Br, I).
    # The '~' operator matches any bond order (single, aromatic etc).
    carbon_halogen_smarts = "[#6]~[#9,17,35,53]"  # #9: F, #17: Cl, #35: Br, #53: I
    halogen_pattern = Chem.MolFromSmarts(carbon_halogen_smarts)
    
    # Check for the carbon-halogen bond in the molecule.
    if mol.HasSubstructMatch(halogen_pattern):
        return True, "Molecule contains at least one carbon-halogen bond"
    else:
        return False, "No carbon-halogen bond found in the molecule"

# Example usage and testing (can be removed/commented out in production)
if __name__ == "__main__":
    test_smiles = [
        "COc1ccc(cc1)C(=O)C(Br)CS(=O)(=O)c1ccc(C)cc1",   # has a Br (bromine) bond
        "Oc1c(Cl)cc(Cl)cc1-c1ccccc1",                      # has Cl (chlorine) bonds
        "CCO",                                            # no halogen bonds
    ]
    
    for smi in test_smiles:
        result, reason = is_organohalogen_compound(smi)
        print(f"SMILES: {smi}\nResult: {result}, Reason: {reason}\n")