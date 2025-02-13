"""
Classifies: CHEBI:24043 flavones
"""
"""
Classifies: Flavones (compounds with a 2-aryl-1-benzopyran-4-one scaffold and substituted derivatives)
We try to require that the substructure match covers the entire core.
"""

from rdkit import Chem

def is_flavones(smiles: str):
    """
    Determines if a molecule is a flavone (or substituted derivative) based on its SMILES string.
    A flavone is defined by the presence of a 2-aryl-1-benzopyran-4-one core.
    
    This improved version uses a SMARTS pattern to identify the core, and then verifies
    that the match covers exactly the number of atoms expected for the core.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule contains the complete flavone core, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern for the flavone core.
    # This pattern encodes a 2-arylchromen-4-one (flavone) skeleton:
    # - "c1ccc(cc1)" represents the 2-aryl (phenyl) group;
    # - "c2oc3ccccc3c(=O)c2" represents the benzopyran-4-one fused skeleton.
    #
    # NOTE: In the previous iteration, molecules were (mis)classified if a partial match was found.
    # Here, after obtaining all substructure matches we require that one match covers exactly the number
    # of atoms specified by the query pattern.
    flavone_smarts = "c1ccc(cc1)-c2oc3ccccc3c(=O)c2"
    query = Chem.MolFromSmarts(flavone_smarts)
    if query is None:
        return False, "Error creating SMARTS pattern for flavone core"
    
    # Find all substructure matches for the flavone core.
    matches = mol.GetSubstructMatches(query)
    if not matches:
        return False, "Flavone core (2-aryl-1-benzopyran-4-one) not found in the molecule"
    
    # We require that one of the matches covers exactly the number of atoms of the query.
    # This avoids accepting partial matches from extra decorated molecules.
    expected_n_atoms = query.GetNumAtoms()
    for match in matches:
        if len(match) == expected_n_atoms:
            return True, "Molecule contains the flavone core (2-aryl-1-benzopyran-4-one) skeleton"
            
    return False, "Flavone-like fragment found but does not match the full expected core pattern"

# Example usage (for testing purposes):
if __name__ == "__main__":
    # Example: diosmetin (a known flavone) from the provided examples.
    test_smiles = "COC1=CC=C(C=C1O)C1=CC(=O)C2=C(O)C=C(O)C=C2O1"
    result, reason = is_flavones(test_smiles)
    print(f"Result: {result}\nReason: {reason}")