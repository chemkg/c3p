"""
Classifies: CHEBI:25676 oligopeptide
"""
"""
Classifies: A peptide containing a relatively small number of amino acids (oligopeptide)
Based on a simple count of peptide (amide) bonds.
"""

from rdkit import Chem

def is_oligopeptide(smiles: str):
    """
    Determines if a molecule is an oligopeptide based on its SMILES string.
    A peptide (with at least 2 amino acids) having a relatively small number of amino acids (chosen here as <= 10) 
    is considered an oligopeptide.
    
    The determination is based on:
      1. Parsing the SMILES.
      2. Searching for peptide (amide) bonds, defined here as the substructure "C(=O)N".
         In a linear peptide the number of amino acids is (number of amide bonds + 1).
         
    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is considered an oligopeptide, False otherwise.
        str: Reason for the classification.
    """
    
    # Parse the SMILES into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern for a peptide (amide) bond.
    # This simplistic pattern "C(=O)N" may also capture other amides,
    # but in the context of known peptides it is sufficient.
    peptide_bond_smarts = "C(=O)N"
    peptide_bond = Chem.MolFromSmarts(peptide_bond_smarts)
    
    # Find all matches for the peptide bond in the molecule
    matches = mol.GetSubstructMatches(peptide_bond)
    num_peptide_bonds = len(matches)
    
    # In a linear peptide, number of amino acids is num_peptide_bonds + 1.
    # (A cyclic peptide would have equal number of bonds and residues, but for simplicity we assume linear chains.)
    num_residues = num_peptide_bonds + 1  if num_peptide_bonds > 0 else 0

    if num_residues == 0:
        return False, "No peptide (amide) bonds detected; not a peptide."
    if num_residues == 1:
        return False, "Only a single amino acid detected; not a peptide."
    
    # Set a threshold for what we consider a 'small number' of amino acids.
    # Here we choose peptides with 2 to 10 amino acids as oligopeptides.
    if num_residues > 10:
        return False, f"Peptide has {num_residues} amino acids, which is too many to be considered an oligopeptide."
    
    return True, f"Peptide detected with {num_residues} amino acids; classified as an oligopeptide."

# Example test (uncomment the following lines to test with one of the provided SMILES)
# test_smiles = "O=C(N[C@@H](CC1=CC=CC=C1)C(O)=O)[C@@H](NC(=O)[C@@H](N)CCCN=C(N)N)CCCN=C(N)N"  # Arg-Arg-Phe example
# print(is_oligopeptide(test_smiles))