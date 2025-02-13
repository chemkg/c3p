"""
Classifies: CHEBI:87657 octanoate ester
"""
"""
Classifies: Octanoate ester 
Definition: Any fatty acid ester in which the carboxylic acid component is octanoic acid (caprylic acid).
"""
from rdkit import Chem

def is_octanoate_ester(smiles: str):
    """
    Determines if a molecule is an octanoate ester based on its SMILES string.
    An octanoate ester is any fatty acid ester with the acyl group being octanoate (from octanoic acid, caprylic acid).
    
    For octanoic acid, the acyl fragment appears as: CCCCCCCC(=O)O in SMILES.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule contains an octanoate ester moiety, False otherwise.
        str: Reason for the classification.
    """
    # Convert SMILES to a molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS pattern for an octanoate ester acyl group.
    # This pattern represents an 8-carbon chain (octyl, where the terminal carbon is the carbonyl carbon)
    # attached with a carbonyl and then an oxygen, as in: CCCCCCCC(=O)O.
    octanoate_pattern = Chem.MolFromSmarts("CCCCCCCC(=O)O")
    if octanoate_pattern is None:
        return False, "Failed to create octanoate ester pattern"
    
    # Check if the molecule contains at least one octanoate ester substructure.
    if mol.HasSubstructMatch(octanoate_pattern):
        return True, "Molecule contains an octanoate ester moiety"
    else:
        return False, "Molecule does not contain an octanoate ester moiety"

# Example usage:
if __name__ == "__main__":
    test_smiles_list = [
        "CCCCCCCC(=O)OCC",  # ethyl octanoate
        "CCCCCCCC(=O)OC",   # methyl octanoate
        "CC(=O)OCCCCCCCC(=O)O",  # a molecule with one octanoate ester and one acetyl group
        "CCCCCCCC(=O)OC[C@H](O)CO" # 1-octanoyl-sn-glycerol
    ]
    for smi in test_smiles_list:
        result, reason = is_octanoate_ester(smi)
        print(f"SMILES: {smi}\nClassification: {result}\nReason: {reason}\n")