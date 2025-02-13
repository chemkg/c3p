"""
Classifies: CHEBI:35759 1-monoglyceride
"""
"""
Classifies: 1-monoglyceride, a monoglyceride in which the acyl substituent is located at position 1.
A 1-monoglyceride is defined as a glycerol molecule (HO–CH2–CHOH–CH2OH) that has exactly one acyl
group attached as an ester at one of its primary (CH2) positions.
The logic here searches for the specific motif:
    Acyl: C(=O)–O–CH2–CHOH–CH2OH
encoded as the SMARTS: "C(=O)O[CH2]C(O)[CH2]O"
(using useChirality=False so that stereochemical decorations do not interfere).
If exactly one match is found we assume the molecule belongs to the 1-monoglyceride class.
"""

from rdkit import Chem

def is_1_monoglyceride(smiles: str):
    """
    Determines if a molecule (given by its SMILES string) is a 1-monoglyceride.
    A 1-monoglyceride has a glycerol backbone in which one of the primary (CH2) hydroxyl groups 
    is esterified with an acyl chain (the motif: R-C(=O)OCH2CHOHCH2OH).
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a 1-monoglyceride, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES into an RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern for the typical 1-monoglyceride motif:
    # R-C(=O)O[CH2]C(O)[CH2]O
    # where the acyl part, R, is attached to a glycerol unit; note that the motif sidesteps
    # other possible ester functions in the molecule.
    mono_smarts = "C(=O)O[CH2]C(O)[CH2]O"
    mono_pattern = Chem.MolFromSmarts(mono_smarts)
    if mono_pattern is None:
        return False, "Invalid SMARTS pattern"
    
    # Find substructure matches ignoring stereochemistry.
    matches = mol.GetSubstructMatches(mono_pattern, useChirality=False)
    n_matches = len(matches)
    
    if n_matches == 0:
        return False, "1-monoglyceride motif (ester attached to a glycerol CH2) not found"
    elif n_matches > 1:
        return False, f"Multiple potential 1-monoglyceride motifs found ({n_matches} matches)"
    else:
        # Exactly one match found.
        return True, "Found 1-monoglyceride motif: one ester group attached at the primary position of a glycerol backbone"

# Example usage:
if __name__ == "__main__":
    # Test one example: rac-1-monodecanoylglycerol
    test_smiles = "CCCCCCCCCC(=O)OCC(O)CO"
    result, reason = is_1_monoglyceride(test_smiles)
    print(result, reason)