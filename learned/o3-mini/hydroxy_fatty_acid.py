"""
Classifies: CHEBI:24654 hydroxy fatty acid
"""
"""
Classifies: Hydroxy fatty acid
Definition: A fatty acid (carboxylic acid with an aliphatic chain) carrying one or more hydroxy substituents on the chain.
Examples include (6S)-6-hydroxyheptanoic acid, ricinelaidic acid, and many others.
"""

from rdkit import Chem

def is_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is a hydroxy fatty acid based on its SMILES string.
    
    A hydroxy fatty acid is defined as a fatty acid (i.e. a molecule with a carboxylic acid group)
    that also carries one or more hydroxy (OH) substituents which are not part of the carboxylic acid group.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if the molecule is a hydroxy fatty acid, False otherwise.
        str: Explanation of the classification decision.
    """
    
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # Define SMARTS for a carboxylic acid group.
    # The pattern "[CX3](=O)[OX2H1]" finds a carbon (with three connections) that is double bonded to an oxygen (implicitly)
    # and single bonded to an -OH group.
    carboxyl_smarts = "[CX3](=O)[OX2H1]"
    carboxyl_query = Chem.MolFromSmarts(carboxyl_smarts)
    if not mol.HasSubstructMatch(carboxyl_query):
        return False, "No carboxylic acid group found; not a fatty acid."
    
    # Identify the oxygen atoms that are part of a carboxylic acid.
    # Here we extract the matching atoms from the carboxyl group. The SMARTS pattern has two explicit atoms:
    # the carbon and the hydroxyl oxygen. We use the hydroxyl oxygen (at position 1) of each match.
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_query)
    carboxyl_oh_indices = set(match[1] for match in carboxyl_matches)
    
    # Define SMARTS for any hydroxyl group (–OH)
    hydroxyl_smarts = "[OX2H]"
    hydroxyl_query = Chem.MolFromSmarts(hydroxyl_smarts)
    free_hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_query)
    
    # Check if there is any hydroxyl group which is not already part of a carboxylic acid group.
    has_free_hydroxy = False
    for match in free_hydroxyl_matches:
        # match is a tuple with the index of the oxygen atom.
        oh_idx = match[0]
        if oh_idx not in carboxyl_oh_indices:
            has_free_hydroxy = True
            break
    
    if not has_free_hydroxy:
        return False, "Fatty acid found but no additional hydroxy substituent detected."
    
    return True, "Contains a carboxylic acid group and a free hydroxy substituent, classifying it as a hydroxy fatty acid."

# Example usage (uncomment to test):
# print(is_hydroxy_fatty_acid("C[C@H](O)CCCCC(O)=O"))  # (6S)-6-hydroxyheptanoic acid