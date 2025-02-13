"""
Classifies: CHEBI:47908 alkanethiol
"""
"""
Classifies: Alkanethiol
An alkanethiol is defined as a compound in which a sulfanyl group (-SH) 
is attached to an alkyl group (i.e., the sulfur is bonded to at least one carbon).
Examples include ethanethiol (CCS), 1-hexanethiol (SCCCCCC), and many others.
"""
from rdkit import Chem

def is_alkanethiol(smiles: str):
    """
    Determines whether a molecule is an alkanethiol by checking if it contains
    a thiol functional group (-SH) attached to a carbon atom.
    
    The approach is to:
      - Parse the SMILES string.
      - Add explicit hydrogens (this ensures that the –SH is correctly represented).
      - Search for a substructure matching the SMARTS pattern "[#16;X2&H1]-[C]".
      
    Args:
        smiles (str): SMILES string representing the molecule.
    
    Returns:
        bool: True if the molecule contains a valid alkanethiol group, False otherwise.
        str: Reason for the classification decision.
    """
    # Parse the SMILES string to an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens to the molecule
    mol = Chem.AddHs(mol)
    
    # Define a SMARTS pattern for a thiol group: sulfur with one hydrogen connected to any carbon.
    # Explanation:
    #   [#16]       -> any sulfur atom (atomic number 16)
    #   X2           -> with exactly two connections (as typical for R-SH)
    #   &H1          -> with exactly one hydrogen attached
    #   -[C]        -> connected directly to any carbon atom
    thiol_pattern = Chem.MolFromSmarts("[#16;X2&H1]-[C]")
    
    # Check if the molecule has the desired thiol substructure
    if mol.HasSubstructMatch(thiol_pattern):
        return True, "Contains a thiol group (-SH) attached to a carbon (alkyl group)."
    else:
        return False, "No alkanethiol functional group (-SH attached to an alkyl group) found."

# Example test cases (uncomment if you wish to run directly)
# test_smiles = [
#     "Cl.SCCN", "CCS", "SCCCCCCS", "SC(C(C)C=O)CC", "SCCCCCC", "S\\C=C(\\CC)/C"
# ]
# for sm in test_smiles:
#     result, reason = is_alkanethiol(sm)
#     print(f"SMILES: {sm} -> {result}; {reason}")