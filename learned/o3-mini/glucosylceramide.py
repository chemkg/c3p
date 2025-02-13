"""
Classifies: CHEBI:36500 glucosylceramide
"""
"""
Classifies: glucosylceramide 
Definition: Any of the cerebrosides in which the monosaccharide head group is glucose.
The function is_glucosylceramide returns True along with a success reason if the molecule contains
a ceramide backbone (identified by an amide motif) and a beta-D-glucosyl sugar moiety attached;
otherwise, it returns False and a reason for failure.
"""

from rdkit import Chem

def is_glucosylceramide(smiles: str):
    """
    Determines if a molecule is a glucosylceramide based on its SMILES string.
    A glucosylceramide is a cerebroside with a beta-D-glucosyl sugar head group attached to the ceramide.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is a glucosylceramide, False otherwise.
        str: Reason for classification.
    """
    # Parse the SMILES string into an RDKit molecule object.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a simplified SMARTS for the ceramide backbone.
    # Here we look for an amide bond where the nitrogen is attached to a chiral carbon (as is common in ceramide).
    ceramide_pattern = Chem.MolFromSmarts("C(=O)N[C@@H]")
    if not mol.HasSubstructMatch(ceramide_pattern):
        return False, "No ceramide backbone (amide motif) found"
    
    # Define a SMARTS for the beta-D-glucosyl moiety.
    # This is based on the common sugar fragment in the provided examples:
    #   CO[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O
    # Note: This is a simplified representation and might not catch all valid variations.
    glucosyl_pattern = Chem.MolFromSmarts("CO[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O")
    if not mol.HasSubstructMatch(glucosyl_pattern):
        return False, "No beta-D-glucosyl head group found"
    
    # If both motifs are present then assume the molecule is a glucosylceramide.
    return True, "Contains ceramide backbone and beta-D-glucosyl head group"

# Example usage (you may remove or comment these out in production):
if __name__ == "__main__":
    test_smiles = "CCCCCCCCCCCCCCCCCC(=O)N[C@@H](CO[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O)[C@H](O)[C@H](O)CCCCCCCCCCC(C)C"
    result, reason = is_glucosylceramide(test_smiles)
    print(result, reason)