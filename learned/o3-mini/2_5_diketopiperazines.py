"""
Classifies: CHEBI:65061 2,5-diketopiperazines
"""
# 2_5_diketopiperazines.py
"""
Classifies: 2,5-diketopiperazines (Any piperazinone that has a piperazine-2,5-dione skeleton)
"""

from rdkit import Chem

def is_2_5_diketopiperazines(smiles: str):
    """
    Determines if a given molecule (by its SMILES string) contains a piperazine-2,5-dione skeleton.
    This skeleton is defined by a six-membered ring with two amide (C(=O)-N) bonds in positions 2 and 5,
    as seen in cyclic dipeptides.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule has a piperazine-2,5-dione skeleton, False otherwise.
        str: Reason for the classification decision.
    """
    # Parse the input SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # Define a SMARTS pattern corresponding to piperazine-2,5-dione.
    # The simplest representation is "O=C1NC(=O)CN1" which captures the cyclic dipeptide skeleton.
    diketopiperazine_smarts = "O=C1NC(=O)CN1"
    pattern = Chem.MolFromSmarts(diketopiperazine_smarts)
    if pattern is None:
        return False, "Error in creating SMARTS pattern."
    
    # Check if the molecule has the 2,5-diketopiperazine (piperazine-2,5-dione) substructure.
    if mol.HasSubstructMatch(pattern):
        return True, "Molecule contains a piperazine-2,5-dione skeleton."
    else:
        return False, "Molecule does not contain a piperazine-2,5-dione skeleton."
        
# Example usage (you can uncomment the following lines to test with one of the provided SMILES):
# result, reason = is_2_5_diketopiperazines("O=C1CNC(=O)CN1")
# print(result, reason)