"""
Classifies: CHEBI:51689 enone
"""
"""
Classifies: an enone (alpha,beta-unsaturated ketone with R(4) ≠ H)
An enone must have a conjugated C=C—C(=O) moiety.
For example, an enone follows the general formula:
    R(1)R(2)C=CR(3)-C(=O)R(4)
with the requirement that R(4) is not hydrogen.
"""

from rdkit import Chem

def is_enone(smiles: str):
    """
    Determines if a molecule is an enone based on its SMILES string.
    
    An enone is defined as an alpha,beta-unsaturated ketone where the carbonyl is in conjugation
    with a carbon–carbon double bond and has a substituent (R4 ≠ H) at the carbonyl carbon.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if the molecule is an enone, False otherwise
        str: Explanation/reason for the classification
    """
    
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern for the enone substructure.
    # This pattern looks for 4 connected carbon atoms:
    #   [#6]=[#6]-[#6](=O)[#6]
    # which corresponds to: carbon=carbon-singlebond-carbonyl-carbon
    # This ensures that the ketone carbon (C=O) is bonded to a carbon (R(4) is not H).
    enone_smarts = "[#6]=[#6]-[#6](=O)[#6]"
    enone_pattern = Chem.MolFromSmarts(enone_smarts)
    
    # Check if the pattern is present in the molecule (at least one match)
    if mol.HasSubstructMatch(enone_pattern):
        return True, "Enone motif found: alpha,beta-unsaturated ketone with non-hydrogen substituent on the carbonyl"
    else:
        return False, "No enone motif (conjugated C=C–C(=O) with R(4) ≠ H) found in the molecule"

# Example usage (uncomment to test):
# smiles_example = "O=C/C=C/c1ccccc1"  # Example: cinnamaldehyde (an aldehyde, not an enone) -> should be False 
# print(is_enone(smiles_example))