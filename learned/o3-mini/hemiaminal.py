"""
Classifies: CHEBI:73080 hemiaminal
"""
"""
Classifies: Hemiaminal
Definition: “Any organic amino compound that has an amino group and a hydroxy group attached to the same carbon atom.
Hemiaminals are intermediates in the formation of imines by addition of an amine to an aldehyde or ketone; 
those derived from primary amines are particularly unstable.”
"""

from rdkit import Chem

def is_hemiaminal(smiles: str):
    """
    Determines if a molecule is a hemiaminal based on its SMILES string.
    A hemiaminal has a carbon with both a hydroxyl group (-OH) and an amino group (-NH2, -NHR, or -NR2)
    attached to the same carbon atom.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule contains a hemiaminal substructure, False otherwise.
        str: A message describing the result of the classification.
    """
    
    # Parse the SMILES string to create a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string; could not parse molecule"
    
    # Define the substructure pattern for a hemiaminal:
    # A tetrahedral (sp3) carbon ([CX4]) attached to a hydroxyl group ([OX2H]) 
    # and an amino group ([NX3]).
    hemiaminal_pattern = Chem.MolFromSmarts("[CX4]([OX2H])([NX3])")
    
    # Search for the hemiaminal pattern in the molecule
    matches = mol.GetSubstructMatches(hemiaminal_pattern)
    if not matches:
        return False, "No carbon with both hydroxyl and amino substituents found; not a hemiaminal"
    
    # If at least one match is found, return True along with the match details
    return True, f"Found hemiaminal substructure in atoms with indices: {matches}"

# Example usage:
if __name__ == "__main__":
    # Test with an example SMILES for 2-Aminopropanol which is a simple hemiaminal
    test_smiles = "OC(N)CC"
    result, reason = is_hemiaminal(test_smiles)
    print(result, reason)