"""
Classifies: CHEBI:27933 beta-lactam antibiotic
"""
"""
Classifies: Beta-lactam Antibiotic 
Definition: An organonitrogen heterocyclic antibiotic that contains a beta-lactam ring.
A beta-lactam ring is a four-membered cyclic amide containing one nitrogen and a carbonyl group.
"""

from rdkit import Chem

def is_beta_lactam_antibiotic(smiles: str):
    """
    Determines if a molecule is a beta-lactam antibiotic based on its SMILES string.
    The test is done by searching for a beta-lactam ring substructure. This is defined here
    as a four-membered ring containing a carbon with a double-bonded oxygen (carbonyl) and a nitrogen.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a beta-lactam antibiotic, False otherwise.
        str: Reason for the classification.
    
    Note: This is a simplistic approach based solely on the presence of a beta-lactam ring.
    Beta-lactam antibiotics are a diverse class of molecules, and further tests (and expert review)
    may be needed for a confident classification.
    """
    # Parse the input SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Define a SMARTS pattern for a beta-lactam ring.
    # Here we assume a simple four-membered ring with a carbonyl carbon, one nitrogen and two additional carbons.
    # The SMARTS "O=C1[C;r4][N;r4][C;r4]1" looks for a ring (closure number "1") that starts with a carbonyl
    # group (O=C) followed by three ring atoms (a carbon, a nitrogen and a carbon) all asserted to be part of a 4-membered ring (r4).
    beta_lactam_smarts = "O=C1[C;r4][N;r4][C;r4]1"
    beta_lactam_pattern = Chem.MolFromSmarts(beta_lactam_smarts)
    if beta_lactam_pattern is None:
        # (This should not happen, but we check just in case)
        return False, "Error in SMARTS pattern for beta-lactam ring"
        
    # Check if the molecule contains the beta-lactam ring substructure
    if mol.HasSubstructMatch(beta_lactam_pattern):
        return True, "Molecule contains a beta-lactam ring, a key structural feature of beta-lactam antibiotics."
    else:
        return False, "Molecule does not contain a beta-lactam ring."