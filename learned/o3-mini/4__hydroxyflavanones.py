"""
Classifies: CHEBI:140331 4'-hydroxyflavanones
"""
"""
Classifies: 4'-hydroxyflavanones
Definition: Any hydroxyflavanone having a hydroxy substituent located at position 4' (i.e. on the para position of the B ring)
"""

from rdkit import Chem

def is_4__hydroxyflavanones(smiles: str):
    """
    Determines if a molecule is a 4'-hydroxyflavanone based on its SMILES string.
    
    A flavanone is characterized by a 2-phenylchroman-4-one (flavanone) core.
    For the subclass of 4'-hydroxyflavanones, we require that the attached phenyl (B ring)
    has a hydroxy substituent at the para position relative to its attachment.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule qualifies as a 4'-hydroxyflavanone, False otherwise.
        str: Reason for classification.
    """
    # Parse the SMILES string into a molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern for the flavanone core.
    # The pattern "C1CC(=O)c2ccccc2O1" approximates the 2-phenylchroman-4-one skeleton.
    core_pattern = Chem.MolFromSmarts("C1CC(=O)c2ccccc2O1")
    if not mol.HasSubstructMatch(core_pattern):
        return False, "Flavanone core (2-phenylchroman-4-one) not found"
    
    # Define a SMARTS pattern to capture the B ring (phenyl group) with a hydroxy group at the para position.
    # We look for a non-aromatic (sp3) carbon attached to an aromatic ring that, in turn, has an -OH substituent.
    # The pattern "[C;!a]c1ccc(O)cc1" approximates a sp3 carbon connected to a benzene ring 
    # which has an -OH group in the position para to the connection.
    b_ring_pattern = Chem.MolFromSmarts("[C;!a]c1ccc(O)cc1")
    if not mol.HasSubstructMatch(b_ring_pattern):
        return False, "No 4'-hydroxy substituent on the B ring found"
    
    return True, "Contains flavanone core with a hydroxy substituent at the 4' position"