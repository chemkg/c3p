"""
Classifies: CHEBI:38958 indole alkaloid
"""
"""
Classifies: Indole alkaloid
Definition: An alkaloid containing an indole skeleton.
"""
from rdkit import Chem

def is_indole_alkaloid(smiles: str):
    """
    Determines if a molecule is an indole alkaloid based on its SMILES string.
    An indole alkaloid is defined as an alkaloid containing an indole skeleton.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as an indole alkaloid, False otherwise.
        str: Reason for classification.
    """
    # Parse the SMILES string into an RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the SMARTS pattern for the indole skeleton.
    # The pattern "c1ccc2c(c1)[nH]cc2" represents a benzene ring fused to a pyrrole ring.
    indole_smarts = "c1ccc2c(c1)[nH]cc2"
    indole_pattern = Chem.MolFromSmarts(indole_smarts)
    
    # Check for the presence of the indole substructure
    if not mol.HasSubstructMatch(indole_pattern):
        return False, "No indole skeleton found"
    
    # Optionally, verify that the molecule contains at least one nitrogen atom,
    # which is required by the definition of an alkaloid.
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    if n_count == 0:
        return False, "No nitrogen atoms found; not an alkaloid"
    
    # Passed both the indole skeleton and nitrogen criteria
    return True, "Molecule contains an indole skeleton and nitrogen, consistent with an indole alkaloid"