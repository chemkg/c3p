"""
Classifies: CHEBI:23824 diol
"""
#!/usr/bin/env python3
"""
Classifies compounds as diols.
A diol is defined as a compound that contains exactly two hydroxy groups.
Examples of diols include (R,R)-butane-2,3-diol, 1,8-tetradecanediol, etc.
"""

from rdkit import Chem

def is_diol(smiles: str):
    """
    Determines if a molecule is a diol based on its SMILES string.
    A diol contains exactly 2 hydroxy groups (–OH).

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is a diol, False otherwise.
        str: Reason for classification.
    """
    # Parse the SMILES string into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a SMARTS pattern for a hydroxy group.
    # [OX2H] represents an oxygen atom (with 2 connections) with an attached hydrogen.
    hydroxy_pattern = Chem.MolFromSmarts("[OX2H]")
    if hydroxy_pattern is None:
        return False, "Failed to create hydroxy group pattern"
    
    # Find all substructure matches to the hydroxy group
    matches = mol.GetSubstructMatches(hydroxy_pattern)
    n_hydroxy = len(matches)

    # Check if we have exactly two hydroxy groups
    if n_hydroxy == 2:
        return True, "Contains exactly two hydroxy groups (diol)"
    else:
        return False, f"Found {n_hydroxy} hydroxy group(s); a diol requires exactly 2"

# Example usage (you can remove this part if only the function is needed)
if __name__ == "__main__":
    example_smiles = "C[C@@H](O)[C@@H](C)O"  # (R,R)-butane-2,3-diol
    result, reason = is_diol(example_smiles)
    print("Is diol?:", result)
    print("Reason:", reason)