"""
Classifies: CHEBI:59644 oxo fatty acid
"""
"""
Classifies: Oxo Fatty Acid
Definition: Any fatty acid containing at least one aldehydic or ketonic group in addition to 
            the carboxylic acid group. Also, the molecule should show a long (predominantly acyclic)
            carbon chain typical of fatty acids.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_oxo_fatty_acid(smiles: str):
    """
    Determines if a molecule is an oxo fatty acid based on its SMILES string.
    An oxo fatty acid is defined as a fatty acid (i.e. possessing a carboxylic acid group,
    and predominantly acyclic long carbon chain) that also contains at least one additional 
    aldehyde or ketone group (i.e. an oxo group) that is not part of the acid group.
    
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
    
    # -------------------------------
    # Check for carboxylic acid group.
    # The SMARTS pattern looks for a carbon with a double-bonded oxygen and a single-bonded hydroxyl.
    acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    acid_matches = mol.GetSubstructMatches(acid_pattern)
    if not acid_matches:
        return False, "No carboxylic acid group detected; not a fatty acid"
    
    # Get the set of carboxyl carbon atom indices (the first atom in the pattern).
    acid_carbon_indices = {match[0] for match in acid_matches}
    
    # -------------------------------
    # Identify additional oxo groups outside the carboxylic acid.
    # SMARTS for ketone: a carbonyl (C=O) with carbons on both sides.
    ketone_pattern = Chem.MolFromSmarts("[#6][CX3](=O)[#6]")
    # SMARTS for aldehyde: a carbonyl (C=O) where the carbon bears a hydrogen.
    aldehyde_pattern = Chem.MolFromSmarts("[#6][CX3H](=O)")
    
    ketone_matches = mol.GetSubstructMatches(ketone_pattern)
    aldehyde_matches = mol.GetSubstructMatches(aldehyde_pattern)
    
    additional_oxo_found = False
    # In both patterns, the carbonyl carbon appears as the second atom (index 1) in the SMARTS.
    for match in ketone_matches:
        if match[1] not in acid_carbon_indices:
            additional_oxo_found = True
            break
    if not additional_oxo_found:
        for match in aldehyde_matches:
            if match[1] not in acid_carbon_indices:
                additional_oxo_found = True
                break
    if not additional_oxo_found:
        return False, "No additional oxo group (aldehyde or ketone) detected outside of the acid"
    
    # -------------------------------
    # Check if the molecule appears "fatty" by verifying the aliphatic (acyclic) carbon chain.
    # We require that the molecule has a decent number of carbon atoms and that the majority (>=70%) 
    # of these carbons are not in a ring (i.e. the molecule is largely acyclic as expected in fatty acids).
    total_carbons = 0
    acyclic_carbons = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:
            total_carbons += 1
            if not atom.IsInRing():
                acyclic_carbons += 1
    if total_carbons < 8:
        return False, "Too few carbon atoms to be considered a fatty acid"
    
    if acyclic_carbons / total_carbons < 0.7:
        return False, "Significant ring content detected; not a typical acyclic fatty acid chain"
    
    # Optionally, one could add more checks (e.g., minimum molecular weight or requirement
    # for a long contiguous alkyl chain) to further reduce ambiguous cases.
    
    return True, "Contains a carboxylic acid group, an additional oxo group outside that acid, and a predominantly acyclic carbon chain typical of a fatty acid"

# Example usage:
if __name__ == "__main__":
    # Using one example: (6E,8Z)-5-oxooctadecadienoic acid
    test_smiles = "O=C(CCCC(O)=O)/C=C/C=C\\CCCCCCCCC"
    result, reason = is_oxo_fatty_acid(test_smiles)
    print("Test molecule classification:", result)
    print("Reason:", reason)