"""
Classifies: CHEBI:23824 diol
"""
"""
Classifies: diol
A diol, for our purposes, is defined as a compound that contains exactly two alcoholic hydroxy groups.
An alcoholic hydroxy group is an –OH group in which the oxygen is directly bonded to an sp³ carbon that is
not part of a carbonyl double bond (i.e. not attached to a carbon that is double‐bonded to oxygen).
In addition, the two –OH groups should be “close” in the molecular framework (no more than 10 bonds apart).
This function uses a SMARTS pattern that restricts the match to –OH groups bonded to sp³ carbons lacking
a C=O bond.
"""

from rdkit import Chem
from rdkit.Chem import rdmolops

def is_diol(smiles: str):
    """
    Determines if a molecule is a diol based on its SMILES string.
    Our criteria: exactly two alcoholic hydroxy groups (–OH groups whose oxygen is attached to a sp³ carbon
    that is not double-bonded to oxygen), and the carbons bearing these –OH groups are within 10 bonds of each other.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if classified as a diol, False otherwise.
        str: Reason for the classification.
    """
    
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern for alcoholic hydroxy groups.
    # This pattern matches an sp³ carbon (CX4) that is not double‐bonded to oxygen (using the exclusion $([CX4]=[OX1]))
    # bonded to an oxygen which carries one hydrogen ([OX2H]).
    alc_oh_smarts = "[CX4;!$([CX4]=[OX1])]-[OX2H]"
    alc_oh_query = Chem.MolFromSmarts(alc_oh_smarts)
    if alc_oh_query is None:
        return False, "Failed to create alcoholic hydroxyl SMARTS pattern"
    
    # Find substructure matches for alcoholic hydroxyl groups.
    # Each match returns a tuple (c_idx, o_idx): the index of the carbon and its attached -OH oxygen.
    alc_matches = mol.GetSubstructMatches(alc_oh_query)
    count_alc = len(alc_matches)
    
    # For a diol we want exactly two alcoholic hydroxyl groups.
    if count_alc != 2:
        return False, f"Found {count_alc} alcoholic hydroxy group(s); a diol requires exactly 2"
    
    # Extract the indices of the carbons attached to the -OH groups.
    alcohol_carbon_indices = [match[0] for match in alc_matches]
    
    # Check the connectivity: the two carbons should not be too far apart in the molecular graph.
    try:
        path = rdmolops.GetShortestPath(mol, alcohol_carbon_indices[0], alcohol_carbon_indices[1])
    except Exception as e:
        return False, f"Error computing shortest path between alcohol-bearing carbons: {e}"
    
    bond_distance = len(path) - 1
    max_bond_distance = 10
    if bond_distance > max_bond_distance:
        return False, f"Alcoholic hydroxy groups are too far apart (distance {bond_distance} bonds; threshold is {max_bond_distance})"
    
    return True, "Contains exactly two alcoholic hydroxy groups (diol) with appropriate connectivity"

# Example usage (you can remove or update this section if only the function is needed):
if __name__ == "__main__":
    # Test with (R,R)-butane-2,3-diol
    example_smiles = "C[C@@H](O)[C@@H](C)O"
    result, reason = is_diol(example_smiles)
    print("Is diol?:", result)
    print("Reason:", reason)