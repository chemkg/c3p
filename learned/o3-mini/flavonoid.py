"""
Classifies: CHEBI:47916 flavonoid
"""
#!/usr/bin/env python
"""
Classifies: Flavonoid
Definition: A flavonoid is any compound with a 1‐benzopyran core (chromene) 
that bears an aryl substituent (typically a phenyl ring) at position 2.
Note: Many natural (and synthetic) flavonoids occur as either the fully aromatic (flavone/flavonol) 
or the partially saturated (flavanone/flavan) forms. 
This program uses two SMARTS patterns – one for an aromatic benzopyran core and one for the partly‐saturated type.
"""

from rdkit import Chem

def is_flavonoid(smiles: str):
    """
    Determines if a molecule is a flavonoid based on its SMILES string.
    The definition used is that a flavonoid must contain a 1-benzopyran 
    (chromene) core having an aryl (phenyl) substituent at the position next to the ring oxygen.
    
    Two SMARTS patterns are applied:
      1. The aromatic (flavone/flavonol) type: the core is fully aromatic.
         SMARTS: "c1ccc2c(c1)oc(c2)~c3ccccc3"
      2. The dihydro (flavanone/flavan) type: the C2 position is sp3.
         SMARTS: "c1ccc2c(c1)OC(C2)~c3ccccc3"
    The "~" bond symbol is used to allow any bond order between the core and the aryl substituent.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule is classified as a flavonoid, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string to an RDKit mol object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS for a fully aromatic benzopyran with an aryl substituent at position 2.
    # Explanation:
    #   "c1ccc2c(c1)oc(c2)" defines a 1-benzopyran (chromene) core.
    #   The terminal fragment "~c3ccccc3" requires that one of the core atoms (at position 2) is directly bonded
    #   (with any bond order, as allowed by "~") to a phenyl ring.
    pattern_flavonoid_arom = Chem.MolFromSmarts("c1ccc2c(c1)oc(c2)~c3ccccc3")
    
    # Define SMARTS for a dihydro (flavanone/flavan) benzopyran type.
    # Here the oxygen-bearing ring contains one sp3 carbon (C2) at the attachment point.
    pattern_flavanone = Chem.MolFromSmarts("c1ccc2c(c1)OC(C2)~c3ccccc3")
    
    # First, check if the molecule matches the aromatic flavonoid pattern.
    if mol.HasSubstructMatch(pattern_flavonoid_arom):
        return True, "Matches aromatic flavonoid skeleton: 1-benzopyran with phenyl substituent at position 2"
    # Next, check for the dihydro (flavanone/flavan) pattern.
    elif mol.HasSubstructMatch(pattern_flavanone):
        return True, "Matches flavanone/flavan skeleton: 1-benzopyran (partially saturated) with phenyl substituent at position 2"
    else:
        return False, "Does not match common flavonoid skeleton patterns"

# Example usage:
if __name__ == "__main__":
    # Test with one known structure (azaleatin)
    test_smiles = "COc1cc(O)cc2oc(-c3ccc(O)c(O)c3)c(O)c(=O)c12"  # azaleatin
    is_fla, reason = is_flavonoid(test_smiles)
    print("Is flavonoid?", is_fla)
    print("Reason:", reason)