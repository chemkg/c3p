"""
Classifies: CHEBI:15841 polypeptide
"""
from rdkit import Chem

def is_polypeptide(smiles: str):
    """
    Determines if a molecule is a polypeptide based on its SMILES string.
    A polypeptide is defined as a peptide containing ten or more amino acid residues.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a polypeptide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Identify peptide bond pattern: N-C(=O)
    peptide_bond_pattern = Chem.MolFromSmarts("N-C(=O)")
    peptide_matches = mol.GetSubstructMatches(peptide_bond_pattern)

    # Count peptide bonds
    num_peptide_bonds = len(peptide_matches)

    # Determine if the molecule is a polypeptide based on the number of peptide bonds
    if num_peptide_bonds >= 9:
        return True, f"Contains {num_peptide_bonds + 1} amino acid residues, qualifies as a polypeptide"
    else:
        return False, f"Contains only {num_peptide_bonds + 1} amino acid residues, not a polypeptide"

    return None, None