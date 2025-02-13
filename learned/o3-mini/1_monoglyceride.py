"""
Classifies: CHEBI:35759 1-monoglyceride
"""
"""
Classifies: 1-monoglyceride, a monoglyceride in which the acyl substituent is located at position 1.
A 1-monoglyceride is defined by having exactly one ester group in which the fatty acyl chain is attached
to the primary (CH2) hydroxyl of a glycerol backbone (HO–CH2–CHOH–CH2OH, with one OH replaced by an acyl group).
"""
from rdkit import Chem

def is_1_monoglyceride(smiles: str):
    """
    Determines if a molecule, given by its SMILES string, is a 1-monoglyceride.
    A 1-monoglyceride has exactly one ester group in which the acyl chain
    is attached at the primary (CH2) position of the glycerol backbone.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule matches 1-monoglyceride criteria, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES into an RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # First, count the number of ester groups in the molecule.
    # An ester group is represented here by the SMARTS pattern [CX3](=O)O.
    ester_pattern = Chem.MolFromSmarts("[CX3](=O)O")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 1:
        return False, f"Expected exactly 1 ester group, found {len(ester_matches)}"
    
    # Next, attempt to identify the 1-monoglyceride motif.
    # For a 1-monoglyceride the acyl chain is attached to the primary (CH2) position of the glycerol backbone.
    # The glycerol backbone in a 1-monoglyceride will appear as R-C(=O)OCH2CH(OH)CH2OH.
    # We use a SMARTS pattern that encodes:
    # "C(=O)O" for the ester functionality and "O[CH2]C(O)C(O)" for the glycerol part.
    mono_pattern = Chem.MolFromSmarts("C(=O)O[CH2]C(O)C(O)")
    if not mol.HasSubstructMatch(mono_pattern):
        return False, "Glycerol backbone with ester at a primary (CH2) position not found"
    
    # If both tests pass then the molecule fits our definition.
    return True, "Found 1-monoglyceride motif: one ester group attached at the primary position of a glycerol backbone"
    
# Example usage:
if __name__ == "__main__":
    # Test one of the provided examples, e.g., rac-1-monodecanoylglycerol.
    test_smiles = "CCCCCCCCCC(=O)OCC(O)CO"
    result, reason = is_1_monoglyceride(test_smiles)
    print(result, reason)