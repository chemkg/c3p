"""
Classifies: CHEBI:16180 N-acylglycine
"""
"""
Classifies: N-acylglycine
Definition: An N-acyl-amino acid in which the amino acid specified is glycine.
N-acylglycine has the key motif R-C(=O)-N-CH2-C(=O)O.
"""

from rdkit import Chem

def is_N_acylglycine(smiles: str):
    """
    Determines if a molecule is an N-acylglycine based on its SMILES string.
    An N-acylglycine must contain an acyl group connected to the amino group of glycine,
    resulting in the motif: R-C(=O)-N-CH2-C(=O)O.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule contains an N-acylglycine motif, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern for the N-acylglycine moiety.
    # The pattern is:
    #   C(=O)  => the acyl carbonyl (can be aliphatic or aromatic)
    #          N => the amide nitrogen
    #      [CH2] => the glycine α-carbon, enforced to be methylene (no side-chain)
    #     C(=O)[O;H1,-1] => the carboxylic acid group (which may appear as -OH or deprotonated)
    #
    # This pattern will match molecules containing an R-C(=O)-N-CH2-C(=O)O fragment.
    acylglycine_smarts = "C(=O)N[CH2]C(=O)[O;H1,-1]"
    pattern = Chem.MolFromSmarts(acylglycine_smarts)
    if pattern is None:
        return False, "Error creating SMARTS pattern"
    
    # Check if the molecule contains the N-acylglycine substructure.
    if not mol.HasSubstructMatch(pattern):
        return False, "The molecule does not contain the N-acylglycine moiety (R-C(=O)-N-CH2-C(=O)O)."
    
    # If the pattern is found, we consider the molecule to be an N-acylglycine.
    return True, "Molecule contains the N-acylglycine moiety (R-C(=O)-N-CH2-C(=O)O)."

# Example testing (can be removed or commented out if desired)
if __name__ == "__main__":
    # Test using N-benzoylglycine as an example:
    test_smiles = "OC(=O)CNC(=O)c1ccccc1"
    result, reason = is_N_acylglycine(test_smiles)
    print(f"SMILES: {test_smiles}\nClassification: {result}\nReason: {reason}\n")