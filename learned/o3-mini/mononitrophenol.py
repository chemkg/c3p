"""
Classifies: CHEBI:39362 mononitrophenol
"""
#!/usr/bin/env python3
"""
Classifies: mononitrophenol
Definition: A nitrophenol is a phenol carrying a single nitro substituent at an unspecified position.
The molecule must contain:
  - Exactly one nitro group (SMARTS: "[N+](=O)[O-]")
  - At least one aromatic hydroxyl group (SMARTS: "[c;R][OX2H]")
Additionally, we verify that the nitro group is attached to an aromatic carbon.
"""

from rdkit import Chem

def is_mononitrophenol(smiles: str):
    """
    Determines if a molecule is a mononitrophenol based on its SMILES string.
    
    A mononitrophenol is defined as a molecule that contains a phenol moiety (i.e. a benzene ring with an attached OH group)
    carrying exactly one nitro substituent at an unspecified position. Additionally, the nitro group should be directly attached
    to an aromatic carbon.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a mononitrophenol, False otherwise.
        str: Explanation of the classification decision.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define substructure patterns
    nitro_pattern = Chem.MolFromSmarts("[N+](=O)[O-]")
    phenol_pattern = Chem.MolFromSmarts("[c;R][OX2H]")  # Aromatic carbon directly bonded to an -OH

    # Check the nitro group: there should be exactly one nitro group
    nitro_matches = mol.GetSubstructMatches(nitro_pattern)
    if len(nitro_matches) == 0:
        return False, "No nitro group found"
    if len(nitro_matches) > 1:
        return False, f"Multiple nitro groups found ({len(nitro_matches)}), expected exactly one"
    
    # Verify that the nitro group is attached to an aromatic carbon.
    # The SMARTS match returns a tuple of indices for the nitro group pattern.
    # In our SMARTS, the first atom is the nitrogen.
    nitro_indices = nitro_matches[0]
    nitro_N = mol.GetAtomWithIdx(nitro_indices[0])
    
    # Find at least one neighbor of the nitro N that is an aromatic carbon and not part of the nitro group itself.
    attached_to_aromatic = False
    for neighbor in nitro_N.GetNeighbors():
        # Skip if neighbor is within the nitro group match itself.
        if neighbor.GetIdx() in nitro_indices:
            continue
        if neighbor.GetAtomicNum() == 6 and neighbor.GetIsAromatic():
            attached_to_aromatic = True
            break
    if not attached_to_aromatic:
        return False, "Nitro group is not attached to an aromatic carbon"
    
    # Check for the presence of a phenol moiety, i.e. an aromatic OH group.
    phenol_matches = mol.GetSubstructMatches(phenol_pattern)
    if not phenol_matches:
        return False, "No aromatic hydroxyl (phenol) group detected"
    
    return True, "Molecule is a mononitrophenol: contains a phenol moiety with a single nitro substituent attached to an aromatic carbon"

# Example usage:
if __name__ == '__main__':
    test_smiles = "Oc1ccc(cc1)[N+]([O-])=O"  # 4-nitrophenol
    result, reason = is_mononitrophenol(test_smiles)
    print(result, reason)