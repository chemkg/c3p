"""
Classifies: CHEBI:35359 carboxamidine
"""
"""
Classifies: Carboxamidine containing compounds
Definition: Compounds having the structure RC(=NR)NR2.
This functional group (commonly seen as -C(=NH)NH2) is recognized by
a carbon atom double bonded to one nitrogen and single bonded to another nitrogen.
Note: This is a simplified substructure search that may match related amidine groups.
"""
from rdkit import Chem

def is_carboxamidine(smiles: str):
    """
    Determines if a molecule contains a carboxamidine group based on its SMILES string.
    The carboxamidine group is defined as RC(=NR)NR2 (commonly the -C(=NH)NH2 group).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule contains a carboxamidine group, False otherwise.
        str: Reason for classification.
    """
    # Parse the provided SMILES string into a molecule object.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the SMARTS pattern for the carboxamidine moiety.
    # This pattern looks for a sp2-hybridized carbon double‐bonded to a nitrogen ([NX2])
    # and single‐bonded to another trivalent nitrogen ([NX3]). This covers the typical
    # format of RC(=NH)NH2 and related variants.
    carboxamidine_pattern = Chem.MolFromSmarts("[CX3](=[NX2])[NX3]")
    
    # Search for the carboxamidine substructure in the molecule.
    matches = mol.GetSubstructMatches(carboxamidine_pattern)
    if matches:
        return True, "Carboxamidine group detected via substructure match."
    else:
        return False, "No carboxamidine group detected."

# Example usage (for testing):
if __name__ == "__main__":
    # Testing with benzamidine: SMILES: NC(=N)c1ccccc1
    test_smiles = "NC(=N)c1ccccc1"
    result, reason = is_carboxamidine(test_smiles)
    print(f"SMILES: {test_smiles}\nResult: {result}\nReason: {reason}")