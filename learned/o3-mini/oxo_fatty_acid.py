"""
Classifies: CHEBI:59644 oxo fatty acid
"""
"""
Classifies: Oxo Fatty Acid
Definition: Any fatty acid containing at least one aldehydic or ketonic group in addition to the carboxylic acid group.
"""

from rdkit import Chem

def is_oxo_fatty_acid(smiles: str):
    """
    Determines if a molecule is an oxo fatty acid based on its SMILES string.
    An oxo fatty acid is defined as a fatty acid (i.e. with a carboxylic acid group)
    containing at least one additional aldehyde or ketone group (i.e. an oxo group) that is not part of the acid.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is an oxo fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES string 
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS pattern for carboxylic acid group.
    # This pattern looks for a carbon atom double-bonded to oxygen and single-bonded to an -OH group.
    acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    acid_matches = mol.GetSubstructMatches(acid_pattern)
    if not acid_matches:
        return False, "No carboxylic acid group detected; not a fatty acid"
    
    # Note: We'll use the carboxyl carbon index (the first atom in the match) to help exclude it later.
    acid_carbons = {match[0] for match in acid_matches}
    
    # Define SMARTS patterns for ketone and aldehyde groups.
    # For ketone: a carbonyl with carbon atoms on both sides.
    ketone_pattern = Chem.MolFromSmarts("[#6][CX3](=O)[#6]")
    # For aldehyde: a carbonyl where the carbon has a hydrogen; note that formyl (H-C(=O)) in an acid (i.e. formic acid)
    # is still a fatty acid so we need to ensure that the acid carbon is not counted twice.
    aldehyde_pattern = Chem.MolFromSmarts("[#6][CX3H](=O)")
    
    ketone_matches = mol.GetSubstructMatches(ketone_pattern)
    aldehyde_matches = mol.GetSubstructMatches(aldehyde_pattern)
    
    # Combine matches from ketone and aldehyde patterns.
    # We'll use sets of the carbon index (first atom in the SMARTS match) so that we can ensure they are not the acid carbon.
    additional_oxo_carbons = set()
    for match in ketone_matches:
        # match[0] is the carbon in the carbonyl group.
        if match[0] not in acid_carbons:
            additional_oxo_carbons.add(match[0])
            
    for match in aldehyde_matches:
        if match[0] not in acid_carbons:
            additional_oxo_carbons.add(match[0])
    
    if not additional_oxo_carbons:
        return False, "No additional oxo group (aldehyde or ketone) detected outside of the carboxylic acid group"
    
    # If needed, one may also check that the molecule appears fatty (e.g., has a long hydrocarbon chain)
    # For now we assume the presence of a carboxylic acid is sufficient to designate it as a fatty acid.
    
    return True, "Contains a carboxylic acid group and an additional oxo group (aldehyde or ketone) outside the acid"

# Example usage:
if __name__ == "__main__":
    # Example: (6E,8Z)-5-oxooctadecadienoic acid
    test_smiles = "O=C(CCCC(O)=O)/C=C/C=C\\CCCCCCCCC"
    result, reason = is_oxo_fatty_acid(test_smiles)
    print("Test molecule classification:", result)
    print("Reason:", reason)