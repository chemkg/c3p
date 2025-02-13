"""
Classifies: CHEBI:46895 lipopeptide
"""
#!/usr/bin/env python3
"""
Classifies: Lipopeptide
Definition: A compound consisting of a peptide with an attached lipid.
This version avoids heavy recursion by using SMARTS patterns.
The algorithm checks for:
  (1) at least one peptide (amide) bond (C(=O)N)
  (2) a lipid tail recognized as a contiguous chain of eight carbon atoms
      (allowing for single or double bonds).
"""

from rdkit import Chem

def is_lipopeptide(smiles: str):
    """
    Determines if a molecule is a lipopeptide based on its SMILES string.
    A lipopeptide is defined as a compound that contains a peptide portion
    (indicated by at least one amide bond, C(=O)N) and an attached lipid tail
    (approximated by the presence of at least eight connected carbon atoms).
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a lipopeptide, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # Check for peptide (amide) bonds. Typical peptide bonds have a C(=O)N functionality.
    amide_smarts = "C(=O)N"
    amide_pattern = Chem.MolFromSmarts(amide_smarts)
    peptide_matches = mol.GetSubstructMatches(amide_pattern)
    if not peptide_matches:
        return False, "No amide (peptide) bonds found in the molecule."
    
    # Check for a lipid chain. The heuristic here is to look for a linear chain of at least 8 carbons.
    # We allow for either single or double bonds between carbons.
    # The SMARTS pattern below looks for 8 carbon atoms connected by any bond (the '~' operator).
    lipid_chain_smarts = "[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]"
    lipid_pattern = Chem.MolFromSmarts(lipid_chain_smarts)
    lipid_matches = mol.GetSubstructMatches(lipid_pattern)
    if not lipid_matches:
        return False, "No contiguous chain of at least eight carbon atoms (lipid tail) found."
    
    # If both peptide bonds and a lipid tail are present the molecule is classified as a lipopeptide.
    return True, ("Molecule contains peptide bonds and has a lipid tail (at least one chain of eight connected carbons).")

# Example usage (for testing purposes)
if __name__ == '__main__':
    # Here we test using one of the sample lipopeptide structures (surfactin A).
    surfactin_A = "[H][C@@]1(CCCCCCCC(C)C)CC(=O)N[C@@H](CCC(O)=O)C(=O)N[C@@H](CC(C)C)C(=O)N[C@H](CC(C)C)C(=O)N[C@@H](C(C)C)C(=O)N[C@@H](CC(O)=O)C(=O)N[C@H](CC(C)C)C(=O)N[C@@H](CC(C)C)C(=O)O1"
    result, reason = is_lipopeptide(surfactin_A)
    print("Surfactin A classification:", result, reason)