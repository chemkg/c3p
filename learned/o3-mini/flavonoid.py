"""
Classifies: CHEBI:47916 flavonoid
"""
#!/usr/bin/env python
"""
Classifies: Flavonoid
Definition: A flavonoid is any compound whose skeleton is based on 1‐benzopyran (chromene)
with an aryl substituent at position 2.
This module uses improved SMARTS patterns that anchor the 2-phenyl substitution using recursive SMARTS 
and allow for extra substituents on the rings.
"""

from rdkit import Chem

def is_flavonoid(smiles: str):
    """
    Determines if a molecule is a flavonoid based on its SMILES string.
    A flavonoid is defined as a compound with a 1-benzopyran (chromene) core having an aryl (phenyl) 
    substituent attached at the C2 position.
    
    Two SMARTS patterns are used:
      1. For the fully aromatic (flavone/flavonol-like) scaffold:
         "c1ccc2c(c1)oc(c2)[$(c3ccccc3)]"
         - This pattern requires a benzopyran core (aromatic ring fused to an oxacycle)
           with a substituent that is recognized recursively as a benzene ring.
      2. For the dihydro (flavanone/flavan-like) scaffold:
         "c1ccc2c(c1)OC(C2)[$(c3ccccc3)]"
         - Similar to the aromatic version, but with the C2 carbon saturated.
         
    The recursive SMARTS syntax "$(...)" ensures that the substituent exactly matches a phenyl group,
    while the rest of the pattern allows other substituents to be present.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a flavonoid, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the improved SMARTS for the aromatic flavonoid scaffold:
    # This requires a benzopyran (chromene) system with a phenyl substituent at position 2.
    pattern_arom = Chem.MolFromSmarts("c1ccc2c(c1)oc(c2)[$(c3ccccc3)]")
    
    # Define the improved SMARTS for the dihydro (flavanone/flavan) flavonoid scaffold.
    pattern_dihydro = Chem.MolFromSmarts("c1ccc2c(c1)OC(C2)[$(c3ccccc3)]")
    
    # Check for aromatic flavonoid scaffold match.
    if mol.HasSubstructMatch(pattern_arom):
        return True, "Matches aromatic flavonoid skeleton: 1-benzopyran with phenyl substituent at C2"
    
    # Check for dihydro flavonoid scaffold match.
    elif mol.HasSubstructMatch(pattern_dihydro):
        return True, "Matches dihydro flavonoid skeleton: 1-benzopyran (saturated at C2) with phenyl substituent at C2"
    
    else:
        return False, "Does not match common flavonoid skeleton patterns"

# Example usage (for testing the function):
if __name__ == "__main__":
    # Example flavonoid: azaleatin (should be classified as a flavonoid)
    test_smiles = "COc1cc(O)cc2oc(-c3ccc(O)c(O)c3)c(O)c(=O)c12"  # azaleatin
    is_fla, reason = is_flavonoid(test_smiles)
    print("Is flavonoid?", is_fla)
    print("Reason:", reason)