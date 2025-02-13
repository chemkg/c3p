"""
Classifies: CHEBI:47916 flavonoid
"""
#!/usr/bin/env python
"""
Classifies: Flavonoid
Definition: Any member of the 'superclass' flavonoids whose skeleton is based on 1-benzopyran 
with an aryl substituent at position 2.
"""

from rdkit import Chem

def is_flavonoid(smiles: str):
    """
    Determines if a molecule is a flavonoid based on its SMILES string.
    Flavonoids are defined as compounds having a 1-benzopyran skeleton (chromene)
    with an aryl substituent at position 2.
    
    For this purpose two SMARTS patterns are used:
    1. A "flavan" pattern (saturated at C2–C3) with the required aryl substituent.
    2. A "flavone" pattern (featuring a carbonyl group at position 4) with the aryl group.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule is a flavonoid, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string to a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS patterns for common flavonoid cores.
    # Pattern 1 - flavan type (2-phenylchroman):
    # This pattern describes a fused benzopyran system (chroman) with a phenyl substituent at the carbon adjacent to the oxygen.
    pattern_flavan = Chem.MolFromSmarts("c1ccc2c(c1)O[C](c3ccccc3)c2")
    # Pattern 2 - flavone type (2-phenylchromen-4-one):
    # This pattern is similar to the above but includes a carbonyl group (C(=O)) indicating oxidation at position 4.
    pattern_flavone = Chem.MolFromSmarts("c1ccc2c(c1)oc(=O)c(c2)c3ccccc3")
    
    # Check if the molecule matches either of the flavonoid core patterns.
    if mol.HasSubstructMatch(pattern_flavan):
        return True, "Matches flavan skeleton: 1-benzopyran with aryl substituent at position 2"
    elif mol.HasSubstructMatch(pattern_flavone):
        return True, "Matches flavone skeleton: 1-benzopyran with aryl substituent at position 2 (with a 4-oxo group)"
    else:
        return False, "Does not match common flavonoid skeleton patterns"

# Example usage:
if __name__ == "__main__":
    test_smiles = "COc1cc(O)cc2oc(-c3ccc(O)c(O)c3)c(O)c(=O)c12"  # azaleatin
    is_fla, reason = is_flavonoid(test_smiles)
    print("Is flavonoid?", is_fla)
    print("Reason:", reason)