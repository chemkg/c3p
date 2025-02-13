"""
Classifies: CHEBI:47908 alkanethiol
"""
"""
Classifies: Alkanethiol
An alkanethiol is a compound in which a sulfanyl group (-SH) is attached to an alkyl group.
Examples include ethanethiol (CCS), 1-hexanethiol (SCCCCCC), etc.
"""
from rdkit import Chem

def is_alkanethiol(smiles: str):
    """
    Determines if a molecule is an alkanethiol based on its SMILES string.
    An alkanethiol is defined as a compound with a sulfanyl group (-SH)
    attached to an alkyl group (i.e., the sulfur is bonded to at least one carbon).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is classified as an alkanethiol, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Add explicit hydrogens so we can correctly detect S-H bonds
    mol = Chem.AddHs(mol)

    # Iterate through all atoms in the molecule
    for atom in mol.GetAtoms():
        # Look for a sulfur atom (atomic number 16)
        if atom.GetAtomicNum() == 16:
            # Check if sulfur has at least one hydrogen attached.
            # GetTotalNumHs() includes both implicit and explicit hydrogens.
            if atom.GetTotalNumHs() > 0:
                # Now ensure that the sulfur is connected to at least one carbon atom (atomic number 6)
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetAtomicNum() == 6:
                        return True, "Contains a thiol group (-SH) attached to a carbon (alkyl group)."
    
    return False, "No alkanethiol functional group (-SH attached to an alkyl group) found."