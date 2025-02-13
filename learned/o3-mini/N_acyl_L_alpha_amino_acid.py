"""
Classifies: CHEBI:48927 N-acyl-L-alpha-amino acid
"""
#!/usr/bin/env python3
"""
Classifies: N-acyl-L-alpha-amino acid
Definition: Any L-alpha-amino acid carrying an N-acyl substituent.
An N-acyl-L-alpha-amino acid should have a chiral alpha-carbon attached to a carboxyl group 
and to an acylated amine (i.e. the nitrogen is bonded to a carbonyl group).
"""

from rdkit import Chem

def is_N_acyl_L_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is an N-acyl-L-alpha-amino acid based on its SMILES string.
    
    The function looks for a chiral alpha-carbon that is bonded to an N-acyl group 
    (i.e. with pattern "NC(=O)") and a carboxyl group (either COOH or COO-).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an N-acyl-L-alpha-amino acid, False otherwise.
        str: Reason for classification.
    """
    # Parse SMILES into a molecule structure
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # We define several SMARTS patterns that capture the core structure:
    # A chiral alpha-carbon ([C@H] or [C@@H]) bonded to a nitrogen that forms an amide
    # (i.e. the N is bonded to a carbonyl: NC(=O)[*]) and also bonded to a carboxyl group,
    # represented as C(=O)O (protonated) or C(=O)[O-] (deprotonated).
    pattern1 = Chem.MolFromSmarts("[C@@H](NC(=O)[*])C(=O)O")
    pattern2 = Chem.MolFromSmarts("[C@H](NC(=O)[*])C(=O)O")
    pattern3 = Chem.MolFromSmarts("[C@@H](NC(=O)[*])C(=O)[O-]")
    pattern4 = Chem.MolFromSmarts("[C@H](NC(=O)[*])C(=O)[O-]")

    # Check if any of the patterns match.
    if (mol.HasSubstructMatch(pattern1) or 
        mol.HasSubstructMatch(pattern2) or 
        mol.HasSubstructMatch(pattern3) or 
        mol.HasSubstructMatch(pattern4)):
        return True, "Molecule has a chiral alpha-carbon bonded to an N-acyl group and a carboxyl group."
    else:
        return False, "The pattern for an N-acyl-L-alpha-amino acid was not found in the molecule."

# Example testing (uncomment for local testing)
# test_smiles = "CC(=O)N[C@@H](CC(O)=O)C(O)=O"  # Example: N-acetyl-L-aspartic acid
# result, reason = is_N_acyl_L_alpha_amino_acid(test_smiles)
# print(result, reason)