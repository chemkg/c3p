"""
Classifies: CHEBI:33566 catechols
"""
"""
Classifies: Catechols (compounds containing an o-diphenol component)
A catechol is defined here as any compound that contains an aromatic ring
bearing two hydroxyl groups in adjacent positions (i.e. a 1,2-dihydroxybenzene fragment).
The previous approach used a rigid SMARTS which led to many false classifications.
This implementation instead uses a more general SMARTS pattern:
    "[c]([OX2H])c([OX2H])"
which requires two aromatic carbons (c) connected by an aromatic bond, each bearing an –OH group.
"""

from rdkit import Chem

def is_catechols(smiles: str):
    """
    Determines if a molecule contains an o-diphenol (catechol) component based on its SMILES string.
    The method uses a SMARTS pattern that matches two connected aromatic atoms (c) with an -OH group each.
    
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
    
    # Define an improved SMARTS pattern to capture an o‐dihydroxybenzene motif.
    # The pattern "[c]([OX2H])c([OX2H])" matches two aromatic atoms (c) directly bonded
    # where each has an attached -OH group (represented as [OX2H]). This is the minimal requirement
    # for an o-diphenol (catechol) moiety.
    catechol_smarts = "[c]([OX2H])c([OX2H])"
    catechol_pattern = Chem.MolFromSmarts(catechol_smarts)
    if catechol_pattern is None:
        return None, None  # Should not occur, but safety check.
        
    # Check if the molecule contains at least one match of the catechol substructure.
    if mol.HasSubstructMatch(catechol_pattern):
        return True, "Found adjacent aromatic hydroxyl groups (o-diphenol component)"
    
    # If no match was found, then no catechol motif exists.
    return False, "No o-diphenol substructure found"

# Uncomment the following test block to run some example tests:
# test_smiles = [
#     # True positives (catechols)
#     "OC1=C(O)C=CC=C1CCCC/C=C\\C/C=C\\CCCCCCCC2=C(O)C(O)=CC=C2",  # Gerronemin F
#     "COc1cc(CCc2ccc(O)c(O)c2)cc(O)c1O",  # dendrocandin E
#     "O[C@@H](CC\\C=C\\c1ccccc1)CCc1ccc(O)c(O)c1",  # (-)-(3S)-1-(3,4-dihydroxyphenyl)-7-phenyl-(6E)-6-hepten-3-ol
#     # False example (catechol not present)
#     "C1=CC=CC=C1",  # benzene without -OH groups
# ]
#
# for s in test_smiles:
#     flag, reason = is_catechols(s)
#     print(f"SMILES: {s}\nClassification: {flag}\nReason: {reason}\n")