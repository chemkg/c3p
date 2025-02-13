"""
Classifies: CHEBI:35267 quaternary ammonium ion
"""
"""
Classifies: Quaternary Ammonium Ion
A quaternary ammonium ion is defined as a derivative of ammonium, NH4(+),
in which all four hydrogens have been replaced by univalent (usually organyl) groups.
Thus, a quaternary ammonium group should have a nitrogen atom with a positive formal charge,
four substituents, and no attached hydrogen atoms.
"""
from rdkit import Chem

def is_quaternary_ammonium_ion(smiles: str):
    """
    Determines if a molecule contains a quaternary ammonium ion based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule contains a quaternary ammonium ion, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Iterate over atoms to find a nitrogen that qualifies
    for atom in mol.GetAtoms():
        # Only consider nitrogen atoms
        if atom.GetAtomicNum() == 7:
            # Check formal charge is +1 (typical for quaternary ammonium)
            if atom.GetFormalCharge() == 1:
                # For a quaternary ammonium, nitrogen should have exactly 4 bonds
                if atom.GetDegree() == 4:
                    # Ensure there are no attached hydrogen atoms (all substituents are organyl)
                    # GetTotalNumHs() gives the total number of hydrogens (explicit and implicit)
                    if atom.GetTotalNumHs() == 0:
                        return True, f"Found quaternary ammonium ion at atom index {atom.GetIdx()}"
                    else:
                        return False, f"Nitrogen at atom index {atom.GetIdx()} has attached hydrogens"
    return False, "No quaternary ammonium ion found in the molecule"

# You can run tests by calling the function with example SMILES strings if needed.
if __name__ == "__main__":
    test_smiles = "CC[N+](C)(C)Cc1ccccc1Br"  # bretylium example from the list
    result, reason = is_quaternary_ammonium_ion(test_smiles)
    print(result, reason)