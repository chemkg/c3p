"""
Classifies: CHEBI:15841 polypeptide
"""
"""
Classifies: polypeptide – A peptide containing ten or more amino acid residues.
This version uses a simplified SMARTS for amide bonds to count peptide linkages.
"""

from rdkit import Chem

def is_polypeptide(smiles: str):
    """
    Determines if a molecule is a polypeptide (a peptide with 10 or more amino acid residues)
    based on its SMILES string.
    
    The function adopts a heuristic approach: it looks for amide bonds (using the pattern "N-C(=O)")
    which are typical for peptide backbones. For linear peptides, the residue count is 
    (number of peptide bonds) + 1.
    
    Note: This approach may count side-chain amide functionalities as well. However,
    with long polypeptides the backbone amide count should dominate.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        (bool, str): A tuple with a boolean indicating whether the molecule is a polypeptide,
                     and a message with details on the residue and peptide bond counts.
    """
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern for an amide bond:
    # This pattern looks for a nitrogen atom directly connected to a carbon atom that is double bonded to an oxygen.
    # (We do not filter based on chirality or extra substituents here to keep the match broad.)
    amide_pattern = Chem.MolFromSmarts("N-C(=O)")
    if amide_pattern is None:
        return False, "Error creating SMARTS pattern"
    
    # Find all substructure matches (each match should correspond to one amide bond).
    matches = mol.GetSubstructMatches(amide_pattern)
    peptide_bond_count = len(matches)
    
    # For a linear peptide, each peptide bond links two amino acid residues.
    # Thus, the residue count is peptide_bond_count + 1.
    residue_count = peptide_bond_count + 1
    
    if residue_count >= 10:
        return True, f"Contains {residue_count} amino acid residues (found {peptide_bond_count} peptide bonds)."
    else:
        return False, f"Only {residue_count} amino acid residues detected (found {peptide_bond_count} peptide bonds)."


# Example usage for testing:
if __name__ == "__main__":
    # Test with one of the provided examples (Variochelin A).
    test_smiles = "O=C(N[C@H](C(=O)N[C@H](C(=O)O)CCCN(O)C(=O)C)CCCN(O)C(=O)C)[C@H]1N(C(=O)[C@H](NC(=O)[C@@H](NC(=O)[C@H]([C@H](O)[C@@H](NC(=O)CCCCCCCCCCC)CCCN=C(N)N)C)[C@H](O)C(=O)O)CO)CCC1"
    result, message = is_polypeptide(test_smiles)
    print(result, message)