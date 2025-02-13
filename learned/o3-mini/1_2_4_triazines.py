"""
Classifies: CHEBI:39410 1,2,4-triazines
"""
"""
Classifies: Compounds containing a 1,2,4-triazine skeleton
Definition: Any compound with a 1,2,4-triazine skeleton, in which nitrogen atoms replace carbon 
at positions 1, 2 and 4 of the core benzene ring structure.
Improvement: We use a SMARTS pattern ("n1ncncc1") to capture a six‐membered aromatic ring with the 
required order (N at positions 1, 2, 4 and C at positions 3,5,6). Then, for each match we extract the 
corresponding sub‐molecule and verify (via its ring information and elemental composition) that it is a standalone 
six‐membered ring having exactly 3 nitrogens and 3 carbons. This extra check helps avoid false positives 
from fused or ambiguous systems and also recovers cases (false negatives) in which the triazine is embedded within 
a fused motif.
"""
from rdkit import Chem

def is_1_2_4_triazines(smiles: str):
    """
    Determines if a molecule contains a valid 1,2,4-triazine skeleton.
    The method uses a SMARTS pattern ("n1ncncc1") to capture a six‐membered aromatic ring with nitrogen atoms 
    at positions 1, 2 and 4 and carbon atoms at positions 3,5,6. For each substructure match, we extract the 
    corresponding sub-molecule and verify that:
      - it has exactly 6 atoms,
      - the sub–molecule itself contains a six–membered ring (using its ring info),
      - its elemental composition is exactly three N and three C.
      
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if a valid 1,2,4‐triazine skeleton is found, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for a six-membered aromatic ring with N at positions 1,2,4
    # and C at positions 3,5,6. The pattern "n1ncncc1" enforces the order in the ring.
    pattern = Chem.MolFromSmarts("n1ncncc1")
    if pattern is None:
        return False, "Error in SMARTS pattern definition"

    # Find all substructure matches of the ring pattern in the molecule.
    matches = mol.GetSubstructMatches(pattern)
    if not matches:
        return False, "No 1,2,4-triazine skeleton found"

    # For each match (a tuple of atom indices) we will try to get a sub-molecule corresponding exactly to it.
    for match in matches:
        # Create a sub-molecule from the atoms in the match.
        submol = Chem.PathToSubmol(mol, list(match))
        # We expect the extracted sub-molecule to have exactly 6 atoms.
        if submol.GetNumAtoms() != 6:
            continue

        # Get ring information from the sub-molecule.
        rings_sub = submol.GetRingInfo().AtomRings()
        # Check that the sub-molecule itself has a six-membered ring.
        if not any(len(r) == 6 for r in rings_sub):
            continue

        # Check that the elemental composition is exactly three C's and three N's.
        elems = sorted([atom.GetSymbol() for atom in submol.GetAtoms()])
        if elems != ['C', 'C', 'C', 'N', 'N', 'N']:
            continue

        # If we got here, then the match appears to be a valid (standalone or well–defined) 1,2,4–triazine ring.
        return True, "Contains a 1,2,4-triazine skeleton"

    return False, "No valid, standalone 1,2,4-triazine ring found"