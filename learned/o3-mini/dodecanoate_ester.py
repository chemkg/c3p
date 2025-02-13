"""
Classifies: CHEBI:87659 dodecanoate ester
"""
#!/usr/bin/env python3
"""
Classifies: dodecanoate ester
Defined as: "Any fatty acid ester in which the carboxylic acid component is lauric acid."
Lauric acid is CH3(CH2)10COOH; when esterified the moiety becomes CH3(CH2)10C(=O)O– (attached to an alcohol).
This program checks for the presence of a laurate ester fragment.
"""
from rdkit import Chem

def is_dodecanoate_ester(smiles: str):
    """
    Determines if a molecule (given by its SMILES string) is a dodecanoate ester.
    A dodecanoate ester is a fatty acid ester in which the acyl component is derived from lauric acid.
    Lauric acid (dodecanoic acid) has the structure CH3(CH2)10COOH, and as an ester yields the fragment:
    CH3(CH2)10C(=O)O– bound to an alkyl group.
    
    For classification, we search for the SMARTS pattern "CCCCCCCCCCCC(=O)O[C]" which represents a chain of 12 carbons,
    with a carbonyl and an ester oxygen attached to another carbon (i.e. not an acid group).
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a dodecanoate ester, False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # Define a SMARTS pattern for a laurate ester fragment.
    # "CCCCCCCCCCCC" represents 12 connected carbons.
    # "(=O)O" represents the ester carbonyl and oxygen.
    # "[C]" ensures that the oxygen is connected to a carbon (ester, not free acid).
    laurate_smarts = "CCCCCCCCCCCC(=O)O[C]"
    
    # Create the query molecule
    query = Chem.MolFromSmarts(laurate_smarts)
    if query is None:
        # In the unlikely event that the SMARTS pattern is invalid
        return False, "Error creating SMARTS for laurate ester."
    
    # Check if the molecule contains the laurate ester substructure.
    if mol.HasSubstructMatch(query):
        return True, "Molecule contains a dodecanoate (laurate) ester functionality."
    else:
        return False, "Molecule does not contain a laurate ester moiety."
        
# Example usage:
if __name__ == "__main__":
    # Test with a known dodecanoate ester example:
    test_smiles = "O(CCCCCC(C)C)C(=O)CCCCCCCCCCCC"  # may represent a laurate ester fragment
    result, reason = is_dodecanoate_ester(test_smiles)
    print(f"Result: {result}\nReason: {reason}")