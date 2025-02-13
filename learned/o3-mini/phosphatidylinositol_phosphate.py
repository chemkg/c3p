"""
Classifies: CHEBI:28765 phosphatidylinositol phosphate
"""
#!/usr/bin/env python3
"""
Classifies: phosphatidylinositol phosphate (a phosphoinositide)
Definition: Any member of the phosphoinositide family of compounds,
of which seven occur naturally.
Improved criteria:
  1. The molecule must have at least two ester bonds as a proxy for fatty acid chains.
     We require the pattern "[#6][C](=O)O" to minimize matching free carboxylic acids.
  2. The molecule must contain an inositol ring. We use a slightly relaxed SMARTS pattern
     for a six‐membered aliphatic ring with hydroxyl substituents.
  3. At least one of the inositol ring carbons must have a substituent oxygen that is directly linked to a phosphate group.
     We check that the phosphorus atom has at least one double‐bonded oxygen (i.e. “=O”).
Note: This is still a simplified filter.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_phosphatidylinositol_phosphate(smiles: str):
    """
    Determines if a molecule is a phosphatidylinositol phosphate based on its SMILES string.
    
    A candidate should have:
      1. Two ester bonds representing two fatty acid chains linked via a glycerol backbone.
      2. An inositol ring (a six‐membered aliphatic ring with multiple –OH groups).
      3. A phosphate group directly attached to one of the inositol hydroxyls.
      
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
    
    # 1. Look for at least two ester bonds.
    # Use a pattern that is common in fatty acids linked as esters:
    ester_pattern = Chem.MolFromSmarts("[#6][C](=O)O")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Found less than 2 ester groups (found {len(ester_matches)}); require at least 2 for fatty acid chains"
    
    # 2. Look for an inositol ring.
    # Allow a six-membered aliphatic ring that has multiple hydroxyl substituents.
    # This SMARTS pattern looks for a six-membered ring (R6) with at least 4 oxygen substituents.
    inositol_pattern = Chem.MolFromSmarts("[C;R6]([OX2H])[C;R6]([OX2H])[C;R6]([OX2H])[C;R6]([OX2H])[C;R6][C;R6]")
    inositol_matches = mol.GetSubstructMatches(inositol_pattern)
    if not inositol_matches:
        return False, "No inositol ring (six-membered ring with several OH groups) detected"
    
    # 3. Check that one of the inositol ring carbons has a substituent oxygen attached to a phosphate group.
    phosphate_found = False
    phosphate_pattern = Chem.MolFromSmarts("P(=O)")  # look for phosphorus with a double bond to oxygen
    for match in inositol_matches:
        for atom_idx in match:
            atom = mol.GetAtomWithIdx(atom_idx)
            # Check neighbors of the inositol carbon
            for nbr in atom.GetNeighbors():
                # We want an oxygen substituent not in the ring (i.e. exocyclic).
                if nbr.GetSymbol() == "O" and nbr.GetIdx() not in match:
                    # Check if this oxygen is bonded to a phosphorus that in turn has a P(=O) bond.
                    for nbr2 in nbr.GetNeighbors():
                        if nbr2.GetSymbol() == "P":
                            # Now check if this phosphorus has at least one double bond to oxygen.
                            phosphorus = nbr2
                            dblO_count = 0
                            for p_nbr in phosphorus.GetNeighbors():
                                # Check bond order between P and neighbor:
                                bond = mol.GetBondBetweenAtoms(phosphorus.GetIdx(), p_nbr.GetIdx())
                                if bond and bond.GetBondTypeAsDouble() >= 2 and p_nbr.GetSymbol() == "O":
                                    dblO_count += 1
                            if dblO_count >= 1:
                                phosphate_found = True
                                break
                    if phosphate_found:
                        break
            if phosphate_found:
                break
        if phosphate_found:
            break

    if not phosphate_found:
        return False, "No phosphate group found attached to the inositol ring"
    
    # If all three criteria are met, we conclude this is a phosphatidylinositol phosphate.
    return True, "Contains two acyl chains linked via ester bonds, an inositol ring, and a phosphate attached to the inositol ring"

# Example usage (uncomment for testing):
# test_smiles = "[C@@H]1(C(C(C([C@H](C1O)OP(=O)(O)O)O)O)O)OP(OC[C@](COC(CCCCCCCCCCCCCCCCC)=O)([H])OC(CCCCCCCCCCCCCCC)=O)(O)=O"
# print(is_phosphatidylinositol_phosphate(test_smiles))