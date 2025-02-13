"""
Classifies: CHEBI:36683 organochlorine compound
"""
"""
Classifies: Organochlorine compound
Definition: An organochlorine compound is a compound containing at least one carbon-chlorine bond.
"""

from rdkit import Chem

def is_organochlorine_compound(smiles: str):
    """
    Determines if a molecule is an organochlorine compound based on its SMILES string.
    
    An organochlorine compound must contain at least one carbon-chlorine (C-Cl) bond.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an organochlorine compound, False otherwise.
        str: Explanation of the result.
    """
    # Parse SMILES string into an RDKit molecule object.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a SMARTS pattern to detect a carbon (atomic number 6) bonded to a chlorine (atomic number 17).
    c_cl_pattern = Chem.MolFromSmarts("[#6]-[#17]")
    if c_cl_pattern is None:
        return False, "Error creating SMARTS pattern for carbon-chlorine bond"
    
    # Check for the presence of at least one carbon-chlorine bond using substructure matching.
    if mol.HasSubstructMatch(c_cl_pattern):
        return True, "Molecule contains at least one carbon-chlorine bond."
    else:
        return False, "No carbon-chlorine bonds found."