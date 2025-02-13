"""
Classifies: CHEBI:52639 N-acylsphingosine
"""
"""
Classifies: N-acylsphingosine
Definition: The parent compounds of the ceramide family, composed of sphingosine having an unspecified fatty acyl group attached to the nitrogen.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_N_acylsphingosine(smiles: str):
    """
    Determines if a molecule is an N-acylsphingosine based on its SMILES string.
    An N-acylsphingosine is defined as a sphingosine backbone bearing an acyl group on its amine.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule is an N-acylsphingosine, False otherwise.
        str: Reason for classification.
    """
    # Parse the SMILES string into an RDKit molecule object.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # First requirement: Check for an amide bond (N-acyl linkage).
    # This SMARTS looks for a nitrogen directly bonded to a carbonyl carbon.
    amide_pattern = Chem.MolFromSmarts("NC(=O)")
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "No amide bond (N-acyl linkage) found"
        
    # Second requirement: Identify the sphingosine backbone.
    # A typical sphingosine (present in ceramide parent compounds) has a backbone where
    # an amino group (now acylated) is connected to two adjacent carbons,
    # one bearing a hydroxyl group and one bearing a CH2OH group.
    # Here we use a heuristic SMARTS pattern that looks for a fragment of the form:
    # [C;H1](O)[C;H1](CO)N[C](=O)
    # which represents: a carbon with an -OH group, next a carbon with a -CH2OH group that is bonded to a nitrogen
    # which in turn is bound to a carbonyl group.
    sphingosine_pattern = Chem.MolFromSmarts("[C;H1](O)[C;H1](CO)N[C](=O)")
    if not mol.HasSubstructMatch(sphingosine_pattern):
        return False, "Sphingosine backbone pattern not found"
    
    # If both patterns are present, we assume this is an N-acylsphingosine.
    return True, "Contains a sphingosine backbone with an N-acyl group attached to the nitrogen"

# Example usage:
if __name__ == "__main__":
    # Test with one of the provided examples: N-(omega-hydroxyoctacosanoyl)sphingosine
    test_smiles = "C(CCCCCCCCCC)CC\\C=C\\[C@@H](O)[C@@H](NC(=O)CCCCCCCCCCCCCCCCCCCCCCCCCCCO)CO"
    result, reason = is_N_acylsphingosine(test_smiles)
    print("N-acylsphingosine:" if result else "Not N-acylsphingosine:", reason)