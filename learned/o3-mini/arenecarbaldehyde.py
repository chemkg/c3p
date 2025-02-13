"""
Classifies: CHEBI:33855 arenecarbaldehyde
"""
#!/usr/bin/env python
"""
Classifies: Arenecarbaldehyde
Definition: Any aldehyde in which the carbonyl group is attached to an aromatic moiety.
The function is_arenecarbaldehyde will return True if the molecule contains an aldehyde group
directly attached to an aromatic carbon, otherwise False with a reason.
"""

from rdkit import Chem

def is_arenecarbaldehyde(smiles: str):
    """
    Determines if a molecule is an arenecarbaldehyde based on its SMILES string.
    An arenecarbaldehyde is any aldehyde in which the carbonyl (–C=O) group is attached to an aromatic moiety.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an arenecarbaldehyde, False otherwise.
        str: Reason for classification.
    """
    
    # Parse the SMILES string into an RDKit molecule object.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the SMARTS pattern for an aromatic-aldehyde:
    # Pattern: an aromatic carbon ([c]) attached via a single bond (-) to an aliphatic carbon [CH],
    # which has a double bond (=) to oxygen (O).
    # This captures the definition of an arenecarbaldehyde.
    arenecarbaldehyde_pattern = Chem.MolFromSmarts('[c]-[CH](=O)')
    if arenecarbaldehyde_pattern is None:
        return False, "Error generating SMARTS pattern"

    # Check for a substructure match; if found, the molecule contains an aldehyde attached to an aromatic moiety.
    if mol.HasSubstructMatch(arenecarbaldehyde_pattern):
        return True, "Found an aldehyde group attached to an aromatic moiety."
    else:
        return False, "No aromatic-attached aldehyde group found."

# For testing purposes (uncomment the lines below):
# test_smiles = [
#     "[H]C(=O)c1cccc2ccccc12",     # 1-naphthaldehyde
#     "Cc1cccc(C=O)c1",             # m-tolualdehyde
#     "[O-][N+](=O)c1ccc(C=O)cc1",   # 4-nitrobenzaldehyde
#     "CC(=O)OCc1cocc2c(C=O)ccc12"   # Baldrinal (example ambiguous structure)
# ]
#
# for smi in test_smiles:
#     result, reason = is_arenecarbaldehyde(smi)
#     print(f"SMILES: {smi}\nResult: {result}, Reason: {reason}\n")