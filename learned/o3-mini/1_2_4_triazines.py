"""
Classifies: CHEBI:39410 1,2,4-triazines
"""
"""
Classifies: Compounds containing a 1,2,4-triazine skeleton
Definition: Any compound with a 1,2,4-triazine skeleton, in which nitrogen atoms replace carbon 
at positions 1, 2 and 4 of the core benzene ring structure.
"""

from rdkit import Chem

def is_1_2_4_triazines(smiles: str):
    """
    Determines if a molecule contains a 1,2,4-triazine skeleton.
    A valid 1,2,4-triazine is a six-membered aromatic ring with exactly three nitrogen atoms 
    (at positions 1, 2, and 4) and three carbon atoms. One key feature is that two of the nitrogens 
    are adjacent (positions 1 and 2).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule contains a 1,2,4-triazine skeleton, False otherwise
        str: A reason for the classification
    """

    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a SMARTS pattern for a 1,2,4-triazine:
    # The pattern represents a six-membered aromatic ring with atoms in cyclic order:
    # N (position 1), N (position 2), C (position 3), N (position 4), C (position 5), C (position 6)
    # Using lower-case letters to denote aromatic atoms.
    # This pattern requires an adjacent pair of N's (positions 1 and 2),
    # which is a distinguishing feature of the 1,2,4-triazine over 1,3,5-triazine.
    triazine_pattern = Chem.MolFromSmarts("n1ncncc1")
    if triazine_pattern is None:
        return False, "Failed to compile SMARTS pattern"

    # Check if the molecule contains the 1,2,4-triazine substructure
    if mol.HasSubstructMatch(triazine_pattern):
        return True, "Contains a 1,2,4-triazine skeleton"
    else:
        return False, "No 1,2,4-triazine skeleton found"