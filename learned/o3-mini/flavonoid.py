"""
Classifies: CHEBI:47916 flavonoid
"""
#!/usr/bin/env python
"""
Classifies: Flavonoid
Definition: A flavonoid is any compound whose skeleton is based on 1‐benzopyran (chromene)
with an aryl substituent at position 2.
Improved approach –
  • We define two SMARTS patterns that explicitly embed the aryl (phenyl) substituent 
    at the C2 position of a benzopyran core.
  • One pattern is for the fully aromatic core (flavone/flavonol type):
      "c1ccc2c(c1)oc(c2)c3ccccc3"
  • The other is for the partially saturated core (flavanone/flavan type):
      "c1ccc2c(c1)OC(C2)c3ccccc3"
  • Although these strict SMARTS may miss some flavonoids that have extra substituents, they
    help reduce false positives from related ring systems.
"""

from rdkit import Chem

def is_flavonoid(smiles: str):
    """
    Determines if a molecule is a flavonoid based on its SMILES string.
    A flavonoid is defined here as a compound having a 1-benzopyran (chromene) core 
    with a phenyl (aryl) substituent attached at the C2 position. Two separate SMARTS
    patterns are used:
      1. For the fully aromatic flavonoids (flavone/flavonol): 
         "c1ccc2c(c1)oc(c2)c3ccccc3"
      2. For the dihydro flavonoids (flavanone/flavan): 
         "c1ccc2c(c1)OC(C2)c3ccccc3"
    Note: Because many natural/synthetic flavonoids bear a variety of extra substituents, 
          a balance must be struck between sensitivity and specificity.
    
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

    # Define SMARTS for the fully aromatic flavonoid skeleton:
    # This pattern ensures:
    #   • a benzopyran (chromene) system: a benzene ring fused with an oxacycle.
    #   • a phenyl ring directly attached to the C2 atom.
    pattern_arom_full = Chem.MolFromSmarts("c1ccc2c(c1)oc(c2)c3ccccc3")
    
    # Define SMARTS for the dihydro (flavanone/flavan) flavonoid skeleton:
    # The only difference is that the C2 atom (bearing the substituent) is saturated.
    pattern_dihydro_full = Chem.MolFromSmarts("c1ccc2c(c1)OC(C2)c3ccccc3")
    
    # First, check for an exact match to the aromatic flavonoid scaffold.
    if mol.HasSubstructMatch(pattern_arom_full):
        return True, "Matches aromatic flavonoid skeleton: 1-benzopyran with phenyl substituent at position 2"
    # Next, check for a match to the dihydro flavonoid scaffold.
    elif mol.HasSubstructMatch(pattern_dihydro_full):
        return True, "Matches dihydro flavonoid skeleton: 1-benzopyran (saturated at C2) with phenyl substituent at position 2"
    else:
        return False, "Does not match common flavonoid skeleton patterns"

# Example usage (for testing the function):
if __name__ == "__main__":
    # Test with one known flavonoid structure (azaleatin)
    test_smiles = "COc1cc(O)cc2oc(-c3ccc(O)c(O)c3)c(O)c(=O)c12"  # azaleatin
    is_fla, reason = is_flavonoid(test_smiles)
    print("Is flavonoid?", is_fla)
    print("Reason:", reason)