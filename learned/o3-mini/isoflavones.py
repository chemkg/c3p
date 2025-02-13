"""
Classifies: CHEBI:38757 isoflavones
"""
"""
Classifies: Isoflavones – any isoflavonoid with a 3-aryl-1-benzopyran-4-one (3-aryl-4H-chromen-4-one) skeleton 
and its substituted derivatives.
The algorithm attempts to match a complete isoflavone skeleton as a substructure.
Two alternatives are allowed:
  1. The fully aromatic isoflavone core.
  2. A dihydroisoflavone core (where the chromenone ring is partially saturated).
Both alternatives require that the external aryl (phenyl) group is attached via a single bond.
"""

from rdkit import Chem

def is_isoflavones(smiles: str):
    """
    Determines if a molecule is an isoflavone 
    (i.e. an isoflavonoid with a 3-aryl-1-benzopyran-4-one skeleton or a dihydro derivative)
    based on its SMILES string.
    
    The query uses a combined SMARTS that allows either:
      - The fully aromatic isoflavone core: "c1ccc2c(c1)oc(=O)c(c2)-c3ccccc3"
      - Or its dihydro counterpart: "c1ccc2C(c1)OC(=O)C2-c3ccccc3"
    The "-" connecting the two parts guarantees that the aryl group is linked by a single bond.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is classified as an isoflavone, False otherwise.
        str: Reason for the classification decision.
    """
    # Parse SMILES into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a SMARTS query that matches either the fully aromatic isoflavone core or the dihydro version.
    # The recursive SMARTS syntax $() is used to allow alternatives.
    isoflavone_smarts = ("$({pattern1}),$({pattern2})").format(
        pattern1="c1ccc2c(c1)oc(=O)c(c2)-c3ccccc3",
        pattern2="c1ccc2C(c1)OC(=O)C2-c3ccccc3"
    )
    
    query = Chem.MolFromSmarts(isoflavone_smarts)
    if query is None:
        return False, "Error constructing SMARTS for isoflavone core"

    # Check if the molecule has a matching substructure.
    if mol.HasSubstructMatch(query):
        return True, ("Molecule matches the isoflavone skeleton (3-aryl-1-benzopyran-4-one or dihydro variant) "
                      "with a non-fused, single-bond-attached phenyl substituent")
    
    return False, ("Molecule does not contain the complete isoflavone skeleton "
                   "(3-aryl-1-benzopyran-4-one core or allowable dihydro derivative)")

# Example usage (for testing purposes):
if __name__ == "__main__":
    # Test with a known isoflavone: daidzein (an isoflavone) SMILES.
    test_smiles = "Oc1ccc(cc1)-c1coc2cc(O)ccc2c1=O"
    result, reason = is_isoflavones(test_smiles)
    print("Result:", result)
    print("Reason:", reason)