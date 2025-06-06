"""
Classifies: CHEBI:23763 pyrroline
"""
from rdkit import Chem

def is_pyrroline(smiles: str):
    """
    Determines if a molecule is a pyrroline based on its SMILES string.
    A pyrroline is an organic heteromonocyclic compound with a structure based on a dihydropyrrole 
    (a five-membered ring containing one nitrogen and at least one double bond).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a pyrroline, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define more comprehensive pyrroline pattern: five-membered ring with one nitrogen
    pyrroline_patterns = [
        "C1=CCNC1",  # 1-Pyrroline
        "C1=CNCC1",  # Alternate positions of double bond
        "C1CC=NC1",  # 2-Pyrroline
        "C1C=CCN1",  # 3-Pyrroline
        "C1CCNC1",   # Dihydropyrrole (both double bonds saturated)
        "C1=CN=CC1", # Pyrrolone variations
        "C1N=CCC1",  # Another pyrroline pattern
    ]

    # Check for any pyrroline substructure match
    for pattern in pyrroline_patterns:
        pyrroline_mol = Chem.MolFromSmarts(pattern) # convert SMARTS pattern to molecule
        if pyrroline_mol and mol.HasSubstructMatch(pyrroline_mol):
            return True, f"Contains pyrroline-like structure: match for pattern {pattern}"

    return False, "No pyrroline-like structure detected"