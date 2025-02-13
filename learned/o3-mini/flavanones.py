"""
Classifies: CHEBI:28863 flavanones
"""
"""
Classifies: Flavanones
Definition: Members of the class of flavans with a 3,4-dihydro-2-aryl-2H-1-benzopyran-4-one skeleton and its substituted derivatives.
"""

from rdkit import Chem

def is_flavanones(smiles: str):
    """
    Determines if a molecule is a flavanone (or substituted derivative) based on its SMILES string.
    The molecule must contain the 3,4-dihydro-2-aryl-2H-1-benzopyran-4-one core.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a flavanone, False otherwise.
        str: Reason for the classification.
    """
    
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern for the flavanone core:
    # This pattern represents the 2-phenylchroman-4-one skeleton:
    #   - "O=C1"  : a carbonyl at position 4 of the heterocycle.
    #   - "CC(c2ccccc2)" : a saturated portion with an aromatic substituent at position 2.
    #   - "Oc2ccccc12" : the oxygen is fused to an aromatic ring, completing the benzopyranone.
    flavanone_pattern = Chem.MolFromSmarts("O=C1CC(c2ccccc2)Oc2ccccc12")
    if flavanone_pattern is None:
        return False, "Failed to create the SMARTS pattern for flavanone core"
    
    # Check if the molecule contains the flavanone substructure.
    if not mol.HasSubstructMatch(flavanone_pattern):
        return False, "Core flavanone skeleton (3,4-dihydro-2-aryl-2H-1-benzopyran-4-one) not found"
    
    # If the core is found, we assume the molecule can be classified as a flavanone (or a substituted derivative)
    return True, "Molecule contains the flavanone core skeleton (3,4-dihydro-2-aryl-2H-1-benzopyran-4-one)"