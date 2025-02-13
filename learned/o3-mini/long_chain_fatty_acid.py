"""
Classifies: CHEBI:15904 long-chain fatty acid
"""
"""
Classifies: Long-chain fatty acid (chain length ranging from C13 to C22)

A long-chain fatty acid is defined here as a molecule containing at least one terminal
carboxylic acid group (i.e. –C(=O)O) and whose total number of carbon atoms (as a first approximation
of chain length) lies between 13 and 22, inclusive.
"""
from rdkit import Chem

def is_long_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a long-chain fatty acid based on its SMILES string.
    
    A long-chain fatty acid is defined as a fatty acid with a chain length ranging 
    from C13 to C22. In our implementation, we require that the molecule contains a 
    terminal carboxylic acid group (–C(=O)O) and that the total number of carbon atoms 
    in the molecule falls within 13 and 22. This is a heuristic approximation.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule qualifies as a long-chain fatty acid, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern for a terminal carboxylic acid group:
    # This matches a carbon atom with a double bond O and an -OH group.
    carboxylic_acid = Chem.MolFromSmarts("[CX3](=O)[OX1H]")
    if not mol.HasSubstructMatch(carboxylic_acid):
        return False, "Molecule does not have a terminal carboxylic acid group"
    
    # Count the number of carbon atoms in the molecule
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    
    # Check if the number of carbons falls within the long-chain fatty acid range (C13 to C22)
    if carbon_count < 13:
        return False, f"Carbon count too low: {carbon_count} carbons (< 13 required)"
    if carbon_count > 22:
        return False, f"Carbon count too high: {carbon_count} carbons (> 22 allowed)"
    
    return True, f"Found terminal carboxylic acid group and carbon count of {carbon_count}"

# Examples for testing (uncomment to run):
# examples = [
#     "O(O)[C@H](CCCCC)\\C=C\\CCCCCCCCCC(O)=O",  # Example 1
#     "OC(=O)CCCCCCCCCCCC(O)=O",                 # Possibly tridecanoic acid
#     "OC(=O)CCCCCCCCCCCCCCCCCCCC(O)=O",         # Too long a chain (more than 22 C)
#     "OC(=O)CCC(O)=O",                          # Too short a chain
# ]
# for s in examples:
#     result, reason = is_long_chain_fatty_acid(s)
#     print(f"SMILES: {s}\nResult: {result}, Reason: {reason}\n")