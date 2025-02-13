"""
Classifies: CHEBI:33855 arenecarbaldehyde
"""
#!/usr/bin/env python
"""
Classifies: Arenecarbaldehyde
Definition: Any aldehyde in which the carbonyl group is attached to an aromatic moiety.
This version improves on the previous by explicitly identifying aldehyde carbons (CH=O) and then checking that
the substituent attached (other than the oxygen) is aromatic and part of a ring.
"""

from rdkit import Chem

def is_arenecarbaldehyde(smiles: str):
    """
    Determines if a molecule is an arenecarbaldehyde based on its SMILES string.
    An arenecarbaldehyde is defined as an aldehyde in which the carbonyl (C=O) group is directly attached
    to an aromatic (ring) carbon.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an arenecarbaldehyde, False otherwise.
        str: Explanation for classification decision.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Our strategy:
    # 1. Identify all aldehyde groups via a SMARTS pattern that captures [CX3H1](=O), i.e., a sp2 carbon,
    #    with exactly one hydrogen and a double bond to oxygen.
    aldehyde_smarts = "[CX3H1](=O)"
    aldehyde_pattern = Chem.MolFromSmarts(aldehyde_smarts)
    if aldehyde_pattern is None:
        return False, "Error generating SMARTS pattern for aldehyde"

    matches = mol.GetSubstructMatches(aldehyde_pattern)
    if not matches:
        return False, "No aldehyde group found"

    # 2. For each matching aldehyde group, check that the carbon (first atom in match) has exactly one neighbor
    #    (other than oxygen) and that neighbor is aromatic and belongs to a ring.
    for match in matches:
        # match is a tuple of atom indices; by our SMARTS the first index corresponds to the aldehyde carbon.
        aldehyde_c_idx = match[0]
        aldehyde_atom = mol.GetAtomWithIdx(aldehyde_c_idx)
        # Find neighbors excluding the oxygen.
        aromatic_neighbor_found = False
        for neighbor in aldehyde_atom.GetNeighbors():
            # Skip oxygen atoms (atomic number 8)
            if neighbor.GetAtomicNum() == 8:
                continue
            # Check if the neighbor is aromatic and in a ring.
            if neighbor.GetIsAromatic() and neighbor.IsInRing():
                aromatic_neighbor_found = True
                break
        if aromatic_neighbor_found:
            return True, "Found an aldehyde group where the carbonyl is directly attached to an aromatic ring."
    
    # If no aldehyde group was found having an aromatic neighbor, classify as false.
    return False, "Aldehyde group found, but not directly attached to an aromatic ring."

# Uncomment below lines for quick testing:
# test_smiles = [
#     "[H]C(=O)c1cccc2ccccc12",     # 1-naphthaldehyde, should be True.
#     "Cc1cccc(C=O)c1",             # m-tolualdehyde, should be True.
#     "CN1C2=CC=CC=C2C(=C1SC3=CC=C(C=C3)Cl)C=O",  # a false positive example from previous model.
#     "[O-][N+](=O)c1ccc(C=O)cc1",   # 4-nitrobenzaldehyde, should be True.
# ]
# for smi in test_smiles:
#     result, reason = is_arenecarbaldehyde(smi)
#     print(f"SMILES: {smi}\nResult: {result}, Reason: {reason}\n")