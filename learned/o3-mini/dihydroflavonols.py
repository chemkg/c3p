"""
Classifies: CHEBI:48039 dihydroflavonols
"""
"""
Classifies: Dihydroflavonols
Definition: Any hydroxyflavanone in which a hydroxy group is present at position 3 
of the heterocyclic ring. The required core is a chroman-4-one with an -OH bonded 
to the saturated carbon at position 3.
"""

from rdkit import Chem

def is_dihydroflavonols(smiles: str):
    """
    Determines if a molecule is a dihydroflavonol based on its SMILES string.
    The molecule must contain a chroman-4-one (flavanone) scaffold in which the saturated
    carbon at position 3 carries a free hydroxyl (-OH) group.
    
    This function uses a single SMARTS pattern that defines the core scaffold with 
    explicit atom mapping for positions 2 and 3. The pattern is:
    
      O=C1[C;!a:2][C;!a:3](O)Oc2ccccc12
      
    Here:
      - "O=C1" identifies the carbonyl on ring atom 1.
      - "[C;!a:2]" is the saturated carbon at position 2.
      - "[C;!a:3](O)" is the saturated carbon at position 3 that must be directly bonded 
         to an -OH.
      - "Oc2ccccc12" closes the oxygen-bridged ring fused with an aromatic system.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule matches the dihydroflavonol scaffold; False otherwise.
        str: Explanation for the classification result.
    """
    # Parse SMILES into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS for dihydroflavonol core:
    # We require a chroman-4-one scaffold with an -OH at position 3.
    # Note: The non-aromatic flag (!a) ensures that the two center carbons are saturated.
    core_smarts = "O=C1[C;!a:2][C;!a:3](O)Oc2ccccc12"
    core_pattern = Chem.MolFromSmarts(core_smarts)
    if core_pattern is None:
        return False, "Internal error: could not create SMARTS pattern."

    # Check if the molecule contains the defined substructure
    if mol.HasSubstructMatch(core_pattern):
        return True, "Found chroman-4-one core with a hydroxyl at position 3 (dihydroflavonol)."
    else:
        return False, "Core scaffold or free hydroxyl at position 3 not found."

# For debugging or unit tests, one might include:
# if __name__ == "__main__":
#     test_smiles = "OC1C(Oc2cc(O)cc(O)c2C1=O)c1ccc(O)c(O)c1"  # Example: (-)-taxifolin
#     result, reason = is_dihydroflavonols(test_smiles)
#     print(result, reason)