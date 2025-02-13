"""
Classifies: CHEBI:15841 polypeptide
"""
"""
Classifies: polypeptide - A peptide containing ten or more amino acid residues.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_polypeptide(smiles: str):
    """
    Determines if a molecule is a polypeptide (a peptide with 10 or more amino acid residues)
    based on its SMILES string.
    
    Heuristic:
      - In a peptide, each amide (peptide) bond (C(=O)N) connects two amino acid residues.
      - For a linear peptide, the number of residues equals the number of peptide bonds plus one.
      - Thus, if we can count at least 9 peptide bonds, then the molecule has 10 or more residues.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule is a polypeptide, False otherwise.
        str: Reason for classification.
    """
    # Parse the SMILES string using RDKit
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern for a peptide (amide) bond.
    # This pattern matches a carbonyl group (C=O) directly attached to a nitrogen (N).
    peptide_bond_pattern = Chem.MolFromSmarts("C(=O)N")
    if peptide_bond_pattern is None:
        return False, "Error creating SMARTS pattern"
    
    # Find all substructure matches for the peptide bond pattern
    matches = mol.GetSubstructMatches(peptide_bond_pattern)
    peptide_bond_count = len(matches)
    
    # In a linear peptide, residue_count = peptide_bond_count + 1.
    residue_count = peptide_bond_count + 1
    
    if residue_count >= 10:
        return True, f"Contains {residue_count} amino acid residues (found {peptide_bond_count} peptide bonds)"
    else:
        return False, f"Only {residue_count} amino acid residues detected (found {peptide_bond_count} peptide bonds)"
        
# Example usage:
if __name__ == "__main__":
    # This is a simple test using one of the provided examples (Variochelin A is generally a peptide)
    test_smiles = "O=C(N[C@H](C(=O)N[C@H](C(=O)O)CCCN(O)C(=O)C)CCCN(O)C(=O)C)[C@H]1N(C(=O)[C@H](NC(=O)[C@@H](NC(=O)[C@H]([C@H](O)[C@@H](NC(=O)CCCCCCCCCCC)CCCN=C(N)N)C)[C@H](O)C(=O)O)CO)CCC1"
    result, reason = is_polypeptide(test_smiles)
    print(result, reason)