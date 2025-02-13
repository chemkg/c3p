"""
Classifies: CHEBI:48039 dihydroflavonols
"""
"""
Classifies: Dihydroflavonols
Definition: Any hydroxyflavanone in which a hydroxy group is present at position 3 
of the heterocyclic (chroman-4-one) ring and where position 2 is substituted by an aromatic ring.
Examples include taxifolin, dihydroisorhamnetin, dihydromyricetin, etc.
"""

from rdkit import Chem

def is_dihydroflavonols(smiles: str):
    """
    Determines if a molecule is a dihydroflavonol based on its SMILES string.
    For a positive classification, the molecule must contain the dihydroflavonol core,
      - a chroman-4-one (flavanone) scaffold fused with an aromatic ring (ring A),
      - with the saturated carbon at position 2 bearing an aromatic substituent (B-ring),
      - and with the next saturated carbon (position 3) bearing a free hydroxyl (-OH).
      
    This function uses a substructure SMARTS pattern with explicit atom mappings:
    
      O[C:1]1[C:2](c3ccccc3)[C:3](O)C(=O)C4=CC=CC=C14

    Here:
      • "O[C:1]1" indicates an oxygen bonded to a saturated carbon (atom 1) that starts a ring.
      • "[C:2](c3ccccc3)" requires that the next saturated carbon (atom 2) is directly
         attached to an aromatic ring (the B-ring).
      • "[C:3](O)" ensures that the following saturated carbon (atom 3) carries a hydroxyl group.
      • "C(=O)" represents the carbonyl group at position 4.
      • "C4=CC=CC=C14" completes the second (A) aromatic ring via ring closure.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the dihydroflavonol core is found; False otherwise.
        str: Explanation for the classification result.
    """
    # Parse the SMILES string into an RDKit Mol object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the SMARTS for the dihydroflavonol core.
    # The pattern enforces:
    #   - a fused heterocycle with a carbonyl (the chroman-4-one core),
    #   - a saturated carbon (position 2) that has a direct aromatic substituent, and
    #   - a hydroxyl group on the next saturated carbon (position 3).
    core_smarts = "O[C:1]1[C:2](c3ccccc3)[C:3](O)C(=O)C4=CC=CC=C14"
    core_pattern = Chem.MolFromSmarts(core_smarts)
    if core_pattern is None:
        return False, "Internal error: could not build SMARTS pattern"
    
    # Check if the molecule contains the dihydroflavonol substructure
    if mol.HasSubstructMatch(core_pattern):
        return True, "Found dihydroflavonol core: chroman-4-one scaffold with an aryl substituent at position 2 and a free hydroxyl at position 3."
    else:
        return False, "Core scaffold (chroman-4-one with 2-aryl substitution and free OH at C3) not found."

# (Optionally, the function below might be used for local testing.)
# if __name__ == "__main__":
#     test_smiles = "OC1C(Oc2cc(O)cc(O)c2C1=O)c1ccc(O)c(O)c1"  # Example: (-)-taxifolin
#     result, reason = is_dihydroflavonols(test_smiles)
#     print(result, reason)