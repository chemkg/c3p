"""
Classifies: CHEBI:29256 thiol
"""
"""
Classifies: Thiol
Definition:
An organosulfur compound in which a thiol group, -SH, is attached to a carbon atom of any aliphatic or aromatic moiety.
"""
from rdkit import Chem

def is_thiol(smiles: str):
    """
    Determines if a molecule is a thiol based on its SMILES string.
    A thiol is defined as an organosulfur compound in which a thiol group (-SH)
    is attached to a carbon atom.

    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a thiol, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string to an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Iterate over all atoms to search for a sulfur with at least one hydrogen
    # and bonded to at least one carbon.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 16:  # Sulfur atom
            # Check if the sulfur atom has one or more hydrogen (implicit or explicit)
            if atom.GetTotalNumHs() > 0:
                # Now check that sulfur is bonded to a carbon atom
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetAtomicNum() == 6:  # Carbon atom
                        return True, "Found thiol group: sulfur with attached hydrogen bonded to carbon"
    
    return False, "No thiol group (-SH attached to carbon) found"