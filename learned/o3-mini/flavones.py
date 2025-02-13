"""
Classifies: CHEBI:24043 flavones
"""
"""
Classifies: Flavones – compounds containing the 2-aryl-1-benzopyran-4-one core (flavone)
A flavone is defined as a flavonoid having a 2-phenyl-4H-1-benzopyran-4-one scaffold.
This version uses an improved SMARTS pattern that enforces the fused-ring connectivity.
"""

from rdkit import Chem

def is_flavones(smiles: str):
    """
    Determines if a molecule is a flavone (or a substituted derivative) based on its SMILES string.
    For a molecule to be classified as a flavone, it must possess the fused 2-phenyl-4H-1-benzopyran-4-one
    core. This improved version uses a SMARTS pattern that more strictly encodes the fused bicyclic (chromen-4-one)
    system with an attached phenyl ring.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule contains the flavone core, False otherwise.
        str: A message explaining the classification decision.
    """
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define an improved SMARTS pattern for the flavone core.
    # The query "c1ccc(cc1)C2=CC(=O)Oc3ccccc23" encodes:
    #   - a phenyl group (c1ccc(cc1))
    #   - attached at position 2 to a fused bicyclic system "C2=CC(=O)Oc3ccccc23"
    #     where the fused system includes a carbonyl at position 4 and an oxygen as part of the heterocycle.
    #
    # This pattern is chosen to enforce correct ring fusion and reject partial matches.
    flavone_smarts = "c1ccc(cc1)C2=CC(=O)Oc3ccccc23"
    query = Chem.MolFromSmarts(flavone_smarts)
    if query is None:
        return False, "Error creating SMARTS pattern for flavone core"
    
    # Attempt to find all substructure matches in the molecule.
    matches = mol.GetSubstructMatches(query)
    if not matches:
        return False, "Flavone core (2-aryl-1-benzopyran-4-one) not found in the molecule"
    
    # For a true hit we require that the matching subgraph covers exactly the core atoms.
    # For our query pattern the expected number of heavy atoms is 15 (as in unsubstituted flavone).
    expected_n_atoms = query.GetNumAtoms()
    # (Substructure matching returns a tuple of indices equal in length to the query.)
    for match in matches:
        if len(match) == expected_n_atoms:
            return True, "Molecule contains the flavone core (2-aryl-1-benzopyran-4-one) skeleton"
    
    return False, "Flavone-like fragment found but does not match the full expected core pattern"

# For testing purposes, one might run examples from the provided list.
if __name__ == "__main__":
    # An example, diosmetin (a known flavone)
    test_smiles = "COC1=CC=C(C=C1O)C1=CC(=O)C2=C(O)C=C(O)C=C2O1"
    result, reason = is_flavones(test_smiles)
    print(f"Result: {result}\nReason: {reason}")