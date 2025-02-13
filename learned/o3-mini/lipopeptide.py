"""
Classifies: CHEBI:46895 lipopeptide
"""
#!/usr/bin/env python3
"""
Classifies: Lipopeptide
Definition: A compound consisting of a peptide with attached lipid.
The algorithm checks for at least one amide bond as a marker for the peptide portion,
and for the presence of a long aliphatic chain (8 consecutive carbons) as a simple lipid signature.
"""

from rdkit import Chem

def is_lipopeptide(smiles: str):
    """
    Determines if a molecule is a lipopeptide based on its SMILES string.
    A lipopeptide should contain a peptide portion (amide bonds) and an attached lipid chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is classified as a lipopeptide, False otherwise.
        str: Reason for the classification.
    """
    # Try to parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Find peptide bonds: we use the amide group pattern "C(=O)N"
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    peptide_matches = mol.GetSubstructMatches(amide_pattern)
    if not peptide_matches:
        return False, "No amide (peptide) bonds found"

    # Find long aliphatic chain: look for at least eight consecutive carbon atoms.
    # This is a rough approximation for a lipid chain.
    lipid_pattern = Chem.MolFromSmarts("CCCCCCCC")
    lipid_matches = mol.GetSubstructMatches(lipid_pattern)
    if not lipid_matches:
        return False, "No lipid chain (>=8 consecutive carbons) found"

    # If both features are detected, classify as lipopeptide.
    return True, "Molecule contains a peptide backbone (amide bonds) with attached lipid chain"

# Example usage (for testing purposes)
if __name__ == '__main__':
    # test with one of the provided examples (for instance, surfactin A)
    surfactin_A = "[H][C@@]1(CCCCCCCC(C)C)CC(=O)N[C@@H](CCC(O)=O)C(=O)N[C@@H](CC(C)C)C(=O)N[C@H](CC(C)C)C(=O)N[C@@H](C(C)C)C(=O)N[C@@H](CC(O)=O)C(=O)N[C@H](CC(C)C)C(=O)N[C@@H](CC(C)C)C(=O)O1"
    result, reason = is_lipopeptide(surfactin_A)
    print("Surfactin A classification:", result, reason)