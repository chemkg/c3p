"""
Classifies: CHEBI:33566 catechols
"""
"""
Classifies: Catechols (compounds containing an o‐diphenol component)
A catechol is defined here as any compound that contains an aromatic ring
bearing two hydroxyl groups in adjacent positions (i.e. a 1,2-dihydroxybenzene fragment).
"""

from rdkit import Chem

def is_catechols(smiles: str):
    """
    Determines if a molecule contains a catechol (o-diphenol) component based on its SMILES string.
    The approach uses a SMARTS query that matches a benzene ring with hydroxyl groups on two adjacent carbons.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule contains an o-diphenol component, False otherwise.
        str: Explanation of the classification.
    """
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern for an o-diphenol (catechol) substructure.
    # The pattern "c1c(O)c(O)cc(c1)" matches a benzene ring (six-membered aromatic ring)
    # with two hydroxyl (–OH) substituents on adjacent carbons.
    # (Note: this pattern may not catch every edge case but focuses on the key o-diphenol motif.)
    catechol_smarts = "c1c(O)c(O)cc(c1)"
    catechol_pattern = Chem.MolFromSmarts(catechol_smarts)
    if catechol_pattern is None:
        return None, None  # if the SMARTS pattern fails for any reason
    
    # Check for the presence of the catechol pattern in the molecule.
    if mol.HasSubstructMatch(catechol_pattern):
        return True, "Found aromatic ring with adjacent hydroxyl groups (o-diphenol component)"
    
    # If no match was found, report that no catechol substructure was detected.
    return False, "No o-diphenol substructure found"

# Example tests (uncomment to run)
# test_cases = [
#     # True positives (catechols)
#     "OC1=C(O)C=CC=C1CCCC/C=C\\C/C=C\\CCCCCCCC2=C(O)C(O)=CC=C2",  # Gerronemin F
#     "COc1cc(CCc2ccc(O)c(O)c2)cc(O)c1O",  # dendrocandin E
#     "O[C@@H](CC\\C=C\\c1ccccc1)CCc1ccc(O)c(O)c1",  # (-)-(3S)-1-(3,4-dihydroxyphenyl)-7-phenyl-(6E)-6-hepten-3-ol
#     # False example (expected to not match the catechol SMARTS)
#     "C1=CC=CC=C1",  # benzene without -OH groups
# ]
# for sm in test_cases:
#     flag, reason = is_catechols(sm)
#     print(f"SMILES: {sm}\nClassification: {flag}\nReason: {reason}\n")