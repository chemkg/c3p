"""
Classifies: CHEBI:27325 xanthophyll
"""
"""
Classifies: Xanthophyll
Definition: A subclass of carotenoids consisting of the oxygenated carotenes.
Xanthophylls typically have a long conjugated polyene chain (common to carotenoids) and contain oxygen-based 
functional groups (such as hydroxyls, ketones, or epoxides) that differentiate them from purely hydrocarbon carotenoids.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_xanthophyll(smiles: str):
    """
    Determines if a molecule is a xanthophyll based on its SMILES string.
    
    A xanthophyll is an oxygenated carotenoid and should have:
    - A long conjugated polyene chain (e.g., several conjugated C=C bonds).
    - At least one oxygen atom present.
    - Sufficient number of carbon atoms to show an extended hydrocarbon framework (typical for carotenoids).
    - A relatively high molecular weight.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is likely a xanthophyll, False otherwise.
        str: Reason for classification.
    """
    
    # Step 1: Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Step 2: Check for the presence of oxygen atoms (the molecule must be oxygenated)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 1:
        return False, "No oxygen atoms found; not oxygenated"
    
    # Step 3: Check that the carbon skeleton is large enough (at least 20 carbons as a rough estimate)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 20:
        return False, "Too few carbon atoms to be a carotenoid"
    
    # Step 4: Look for a conjugated polyene chain indicative of carotenoids.
    # Here we use a simple SMARTS pattern for a chain with three consecutive C=C bonds.
    polyene_pattern = Chem.MolFromSmarts("C=C-C=C-C=C")
    if not mol.HasSubstructMatch(polyene_pattern):
        return False, "No sufficiently long conjugated polyene chain found"
    
    # Step 5: Check that the molecular weight is high enough (carotenoids generally are large molecules)
    mw = rdMolDescriptors.CalcExactMolWt(mol)
    if mw < 300:
        return False, "Molecular weight too low for a carotenoid"
    
    return True, "Contains a long polyene chain and oxygen functionalities consistent with xanthophylls"

# Example usage (this line can be removed in production code):
if __name__ == '__main__':
    test_smiles = "CC(\\C=C\\C=C(C)C=C\\C1=C(C)[C@@H](O)CC1(C)C)"  # partial structure example
    result, reason = is_xanthophyll(test_smiles)
    print(result, reason)