"""
Classifies: CHEBI:39410 1,2,4-triazines
"""
"""
Classifies: Compounds containing a 1,2,4-triazine skeleton
Definition: Any compound with a 1,2,4-triazine skeleton, in which nitrogen atoms replace carbon 
at positions 1, 2 and 4 of the core benzene ring structure.
Improvement: We use a SMARTS pattern ("n1ncncc1") to explicitly capture an aromatic six‐membered ring 
with nitrogen atoms at positions 1,2,4 and carbons at positions 3,5,6. To reduce false positives from 
fused or ambiguous systems, we further check if the matching atoms exactly correspond to a six‐membered ring 
from the molecule’s ring information.
"""

from rdkit import Chem

def is_1_2_4_triazines(smiles: str):
    """
    Determines if a molecule contains a valid 1,2,4-triazine skeleton.
    The improved method uses a SMARTS pattern "n1ncncc1" to capture a six-membered aromatic ring
    containing nitrogen atoms at positions 1,2,4 and carbon atoms at positions 3,5,6. In addition,
    we verify that the match corresponds to an actual six-membered ring in the molecule’s ring system.
    
    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule contains a valid 1,2,4-triazine skeleton, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for a 1,2,4-triazine ring:
    # The pattern "n1ncncc1" represents a six-membered aromatic ring with atoms in the following order:
    # position 1: aromatic nitrogen, position 2: aromatic nitrogen, position 3: aromatic carbon,
    # position 4: aromatic nitrogen, position 5: aromatic carbon, position 6: aromatic carbon.
    pattern = Chem.MolFromSmarts("n1ncncc1")
    if pattern is None:
        return False, "Error in SMARTS pattern definition"

    # Find all substructure matches for the SMARTS pattern
    matches = mol.GetSubstructMatches(pattern)
    if not matches:
        return False, "No 1,2,4-triazine skeleton found"

    # Get the six-membered rings from the molecule's ring info (as sets of atom indices)
    rings = [set(r) for r in mol.GetRingInfo().AtomRings() if len(r) == 6]
    
    # Look for a match that exactly corresponds to one of the six-membered rings
    for match in matches:
        if set(match) in rings:
            return True, "Contains a 1,2,4-triazine skeleton"

    return False, "No valid, standalone 1,2,4-triazine ring found"