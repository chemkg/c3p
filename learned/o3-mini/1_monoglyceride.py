"""
Classifies: CHEBI:35759 1-monoglyceride
"""
"""
Classifies: 1-monoglyceride, a monoglyceride in which the acyl substituent is located at a primary (CH2) position.
A 1-monoglyceride is defined as a glycerol molecule (HO–CH2–CHOH–CH2OH) that has exactly one acyl group (R–C(=O))
attached as an ester on one of the primary positions. 

This version first rejects molecules containing phosphorus (which are typical for phospholipids)
and then searches for the specific motif:
    Acyl portion: R–C(=O)
    Glycerol portion: O[CH2]C([OX2H])[CH2][OX2H]

If exactly one match is found we assume the molecule belongs to the 1-monoglyceride class.
"""
from rdkit import Chem

def is_1_monoglyceride(smiles: str):
    """
    Determines if a molecule (given by its SMILES string) is a 1-monoglyceride.
    A 1-monoglyceride has a glycerol backbone (HO–CH2–CHOH–CH2OH) with exactly one acyl
    group attached as an ester at one of its primary (CH2) positions.
    
    Args:
        smiles (str): SMILES string representing the molecule.
    
    Returns:
        bool: True if the molecule is a 1-monoglyceride, False otherwise.
        str: Reason for the classification.
    """
    # Parse SMILES string to a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Filter out molecules that contain a phosphorus atom (common in phospholipids, not in monoacylglycerols)
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 15:  # Phosphorus atomic number
            return False, "Molecule contains phosphorus; likely not a 1-monoglyceride"
    
    # Define a refined SMARTS pattern for the 1-monoglyceride motif.
    # This pattern enforces that the glycerol part carries free hydroxyl groups:
    #   Acyl: C(=O)
    #   Ester bond: O
    #   Glycerol backbone with free OH groups: [CH2]C([OX2H])[CH2][OX2H]
    mono_smarts = "C(=O)O[CH2]C([OX2H])[CH2][OX2H]"
    mono_pattern = Chem.MolFromSmarts(mono_smarts)
    if mono_pattern is None:
        return False, "Invalid SMARTS pattern"
    
    # Find substructure matches ignoring stereochemistry
    matches = mol.GetSubstructMatches(mono_pattern, useChirality=False)
    n_matches = len(matches)
    
    if n_matches == 0:
        return False, "1-monoglyceride motif not found (ester not attached to a glycerol CH2 with free OH groups)"
    elif n_matches > 1:
        return False, f"Multiple potential 1-monoglyceride motifs found ({n_matches} matches) in molecule"
    else:
        return True, "Found 1-monoglyceride motif: one ester group attached at a primary position of a glycerol backbone"

# Example usage:
if __name__ == "__main__":
    # Test one example: rac-1-monodecanoylglycerol
    test_smiles = "CCCCCCCCCC(=O)OCC(O)CO"
    result, reason = is_1_monoglyceride(test_smiles)
    print(result, reason)