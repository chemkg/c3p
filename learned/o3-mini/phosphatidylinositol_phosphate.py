"""
Classifies: CHEBI:28765 phosphatidylinositol phosphate
"""
#!/usr/bin/env python3
"""
Classifies: phosphatidylinositol phosphate (a phosphoinositide)
Definition: Any member of the phosphoinositide family of compounds,
of which seven occur naturally.
This code checks for the presence of:
  - Two acyl chains attached via ester bonds (diacylglycerol).
  - An inositol ring (a non-aromatic six‐membered ring with multiple –OH groups).
  - A phosphate group attached to the inositol ring.
Note: This is a simplified filter.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_phosphatidylinositol_phosphate(smiles: str):
    """
    Determines if a molecule is a phosphatidylinositol phosphate based on its SMILES string.
    
    A candidate should have:
      1. Two ester (acyl chain) groups (the fatty acid tails).
      2. An inositol ring (a six-membered saturated ring with multiple hydroxyl substituents).
      3. A phosphate group attached to the inositol ring.
      
    Args:
        smiles (str): SMILES representation of the molecule.
        
    Returns:
        bool: True if the molecule meets the criteria for a phosphatidylinositol phosphate,
              False otherwise.
        str: Explanation for the classification.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 1. Look for diacylglycerol ester bonds.
    # The ester bond pattern for a fatty acid is approximated as C(=O)O.
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Found less than 2 ester groups (found {len(ester_matches)}); require at least 2 for fatty acid chains"
    
    # 2. Look for an inositol ring.
    # This substructure is a simplified representation of myo-inositol:
    # A six-membered aliphatic ring with hydroxyl groups.
    inositol_pattern = Chem.MolFromSmarts("[C@@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@H](O)1")
    inositol_matches = mol.GetSubstructMatches(inositol_pattern)
    if not inositol_matches:
        return False, "No inositol ring (6-membered ring with multiple OH groups) detected"
    
    # 3. Check that at least one of the inositol rings has a phosphate group attached.
    # We look for an oxygen attached to the inositol ring that in turn is bound to a phosphorus.
    phosphate_attached_to_inositol = False
    for match in inositol_matches:
        # For each atom in the inositol ring (as matched by our pattern)
        for idx in match:
            atom = mol.GetAtomWithIdx(idx)
            for nbr in atom.GetNeighbors():
                if nbr.GetSymbol() == "O":
                    # Check if this oxygen neighbors a phosphorus atom.
                    for nbr2 in nbr.GetNeighbors():
                        if nbr2.GetSymbol() == "P":
                            # Found phosphate (P) attached to an oxygen of the inositol ring.
                            phosphate_attached_to_inositol = True
                            break
                    if phosphate_attached_to_inositol:
                        break
            if phosphate_attached_to_inositol:
                break
        if phosphate_attached_to_inositol:
            break
    
    if not phosphate_attached_to_inositol:
        return False, "No phosphate group found attached to the inositol ring"
    
    # If all three criteria are met, we conclude this is a phosphatidylinositol phosphate.
    return True, "Contains two acyl chains linked by ester bonds, an inositol ring, and a phosphate attached to the inositol ring"
    
# Uncomment below lines for testing with one of the provided example SMILES.
# test_smiles = "[C@@H]1(C(C(C([C@H](C1O)OP(=O)(O)O)O)O)O)OP(OC[C@](COC(CCCCCCCCCCCCCCCCC)=O)([H])OC(CCCCCCCCCCCCCCC)=O)(O)=O"
# print(is_phosphatidylinositol_phosphate(test_smiles))