"""
Classifies: CHEBI:35358 sulfonamide
"""
"""
Classifies: Sulfonamide
Definition: An amide of a sulfonic acid RS(=O)2NR'2
"""

from rdkit import Chem

def is_sulfonamide(smiles: str):
    """
    Determines if a molecule is a sulfonamide based on its SMILES string.
    A sulfonamide is defined as an amide of a sulfonic acid, with the general 
    structure RS(=O)(=O)NR'2, i.e. a sulfonyl group (S with two double bonded O's)
    bonded to a nitrogen.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is classified as a sulfonamide, False otherwise.
        str: Reason for the classification decision.
    """
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for a sulfonamide group.
    # This pattern looks for a sulfur atom (with tetra-coordination) 
    # that is double-bonded to two oxygen atoms and single-bonded to a nitrogen atom.
    sulfonamide_pattern = Chem.MolFromSmarts("[SX4](=O)(=O)[NX3]")
    
    # Check if the sulfonamide functional group is present in the molecule
    if mol.HasSubstructMatch(sulfonamide_pattern):
        return True, "Molecule contains a sulfonamide group: RS(=O)(=O)NR'2"
    else:
        return False, "No sulfonamide group (S(=O)(=O)N) found in the molecule"
    
# Example usage (uncomment to run tests):
# smiles_examples = [
#     "Brc1ccc(c2ccccc12)S(=O)(=O)NCc1ccccn1",  # pyrabactin contains sulfonamide moiety
#     "CC1=CC=CC=C1S(=O)(=O)N2[C@@H]3CN4C(=O)C=CC=C4[C@H]2[C@H]([C@@H]3CO)C(=O)NCCOC"  # Another example
# ]
# for sm in smiles_examples:
#     result, reason = is_sulfonamide(sm)
#     print(f"SMILES: {sm}\nResult: {result}\nReason: {reason}\n")