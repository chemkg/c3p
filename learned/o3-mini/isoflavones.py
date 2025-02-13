"""
Classifies: CHEBI:38757 isoflavones
"""
"""
Classifies: Isoflavones – any isoflavonoid with a 3-aryl-1-benzopyran-4-one skeleton and its substituted derivatives.
"""
from rdkit import Chem

def is_isoflavones(smiles: str):
    """
    Determines if a molecule is an isoflavone (isoflavonoid with a 3-aryl-1-benzopyran-4-one skeleton)
    based on its SMILES string.
    
    An isoflavone should contain the core fused bicyclic system (benzopyran-4-one) and the aromatic (phenyl) 
    substituent at the position corresponding to position 3.

    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is classified as an isoflavone, False otherwise.
        str: Reason for the classification decision.
    """
    # Parse the SMILES string to create an RDKit molecule object.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # SMARTS pattern for the 3-aryl-1-benzopyran-4-one skeleton.
    # Explanation of the pattern:
    # "c1ccc2c(c1)oc(=O)c2c(-c3ccccc3)"
    #  • c1ccc2c(c1)oc(=O)c2   describes the fused bicyclic benzopyran-4-one core,
    #  • c(-c3ccccc3)         demands an aromatic (phenyl) substituent attached (ideally at the 3-position).
    #
    # Note: This pattern accepts substituted derivatives too.
    core_smarts = "c1ccc2c(c1)oc(=O)c2c(-c3ccccc3)"
    core_mol = Chem.MolFromSmarts(core_smarts)
    if core_mol is None:
        return False, "Error in SMARTS pattern construction"

    # Perform substructure matching.
    if mol.HasSubstructMatch(core_mol):
        return True, "Molecule contains the 3-aryl-1-benzopyran-4-one isoflavone skeleton"
    else:
        return False, "Molecule does not contain the 3-aryl-1-benzopyran-4-one isoflavone skeleton"

# Example usage:
if __name__ == "__main__":
    # You can test the function with an example SMILES such as daidzein.
    test_smiles = "Oc1ccc(cc1)-c1coc2cc(O)ccc2c1=O"  # daidzein
    result, reason = is_isoflavones(test_smiles)
    print("Result:", result)
    print("Reason:", reason)