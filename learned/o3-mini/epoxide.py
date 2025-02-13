"""
Classifies: CHEBI:32955 epoxide
"""
"""
Classifies: Epoxide – Any cyclic ether in which the oxygen atom forms part of a 3-membered ring.
"""

from rdkit import Chem

def is_epoxide(smiles: str):
    """
    Determines if a molecule is an epoxide based on its SMILES string.
    An epoxide is defined as any cyclic ether in which the oxygen atom forms part of a 
    3-membered ring (i.e. a three membered ring containing one oxygen and two carbon atoms).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an epoxide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES into a molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the SMARTS pattern for a 3-membered cyclic ether (epoxide)
    # This pattern looks for a 3-membered ring that has one oxygen atom and two carbon atoms.
    epoxide_pattern = Chem.MolFromSmarts('[C;R3][O;R3][C;R3]')
    
    # Use RDKit's substructure search to look for the epoxide pattern
    if mol.HasSubstructMatch(epoxide_pattern):
        return True, "Contains a 3-membered cyclic ether (epoxide) group"
    else:
        return False, "Does not contain a 3-membered cyclic ether (epoxide) group"