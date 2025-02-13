"""
Classifies: CHEBI:17389 2-monoglyceride
"""
"""
Classifies: 2-monoglyceride, defined as 'A monoglyceride in which the acyl substituent is located at position 2.'
"""

from rdkit import Chem

def is_2_monoglyceride(smiles: str):
    """
    Determines if a molecule is a 2-monoglyceride based on its SMILES string.
    A 2-monoglyceride is a glycerol molecule in which the acyl group is attached
    to the secondary (2) position via an ester bond. The expected substructure is
    a glycerol backbone of the form HOCH2-CHOH-CH2OH where the central hydroxyl is
    replaced with an acyl ester (i.e. C(=O)O- attached at the 2-position).
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule is identified as a 2-monoglyceride, False otherwise.
        str: Explanation for classification.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern for a 2-monoglyceride substructure.
    # This pattern looks for an ester group (C(=O)O) attached to a central glycerol carbon.
    # The central carbon is linked to two CH2O (primary alcohol) groups.
    # SMARTS explanation:
    #   [CX3](=O)       --> an sp2 carbon with a carbonyl group (acyl carbon)
    #   O                --> an oxygen linking the acyl group to glycerol
    #   [CH]            --> the central glycerol carbon (should be CH, since its -OH is replaced by ester)
    #   ([CH2]O)        --> one terminal CH2OH group
    #   [CH2]O         --> the other terminal CH2OH group
    ester_glycerol_smarts = "[CX3](=O)O[CH]([CH2]O)[CH2]O"
    pattern = Chem.MolFromSmarts(ester_glycerol_smarts)
    if pattern is None:
        # In case the SMARTS pattern is not compiled, return error.
        return False, "SMARTS pattern could not be compiled"
    
    # Find all matching substructures in the molecule.
    matches = mol.GetSubstructMatches(pattern)
    
    # For a 2-monoglyceride we expect exactly one occurrence of the pattern.
    if len(matches) == 0:
        return False, "2-monoglyceride substructure not found"
    elif len(matches) > 1:
        return False, "Multiple 2-monoglyceride-like substructures found; molecule may be di- or tri-acylated"
    
    # If exactly one match is found, we classify the molecule as a 2-monoglyceride.
    return True, "Molecule contains a glycerol backbone with a single acyl ester at the 2-position"