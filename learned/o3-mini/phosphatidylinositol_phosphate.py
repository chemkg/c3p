"""
Classifies: CHEBI:28765 phosphatidylinositol phosphate
"""
#!/usr/bin/env python3
"""
Classifies: phosphatidylinositol phosphate (a phosphoinositide)
Definition: Any member of the phosphoinositide family of compounds, of which seven occur naturally.
Improved criteria:
  1. The molecule must have at least two ester bonds as a proxy for fatty acid chains.
     We require the pattern "[#6][C](=O)O" to minimize matching free carboxylic acids.
  2. The molecule must contain an inositol ring. We search for a six‐membered ring made entirely of carbons 
     with at least three exocyclic oxygen substituents.
  3. At least one of the oxygens attached to the inositol ring must be linked to a phosphate group,
     i.e. to a phosphorus atom that has at least one double‐bonded oxygen.
Note: This is still a simplified filter.
"""
from rdkit import Chem

def is_phosphatidylinositol_phosphate(smiles: str):
    """
    Determines if a molecule is a phosphatidylinositol phosphate based on its SMILES string.
    
    A candidate should have:
      1. Two ester bonds representing two fatty acid chains linked via a glycerol backbone.
      2. An inositol ring (a six‐membered aliphatic ring made of carbons with several exocyclic -OH or phosphate substituents).
      3. At least one hydroxyl (or phosphate‐substituted) on the inositol ring is further bound to a phosphate group.
      
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
    
    # 1. Check for at least two ester bonds (fatty acid chains linked to a glycerol backbone)
    ester_pattern = Chem.MolFromSmarts("[#6][C](=O)O")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Found less than 2 ester groups (found {len(ester_matches)}); require at least 2 for fatty acid chains"

    # 2. Search for an inositol ring.
    # Look for six-membered rings (non-aromatic) made entirely of carbon atoms.
    ring_info = mol.GetRingInfo()
    inositol_candidate = None  # will store a tuple: (ring_atom_indices, list of (ring atom index, exocyclic oxygen index))
    for ring in ring_info.AtomRings():
        if len(ring) == 6:
            # Check if all atoms in the ring are carbon (atomic num = 6)
            if all(mol.GetAtomWithIdx(idx).GetAtomicNum() == 6 for idx in ring):
                # Now count the exocyclic oxygen substituents on the ring atoms.
                oxy_subs = []
                for idx in ring:
                    atom = mol.GetAtomWithIdx(idx)
                    for nbr in atom.GetNeighbors():
                        if nbr.GetIdx() not in ring and nbr.GetSymbol() == "O":
                            oxy_subs.append( (idx, nbr.GetIdx()) )
                # You can adjust the threshold (e.g. require at least 3 substituents)
                if len(oxy_subs) >= 3:
                    inositol_candidate = (ring, oxy_subs)
                    break
    if not inositol_candidate:
        return False, "No inositol ring (a six-membered carbon ring with at least 3 oxygen substituents) detected"

    # 3. Verify that at least one of the exocyclic oxygens attached to the inositol ring is linked to a phosphate group.
    phosphate_found = False
    # For each oxygen substituent on the inositol ring:
    for ring_atom_idx, oxy_idx in inositol_candidate[1]:
        oxy_atom = mol.GetAtomWithIdx(oxy_idx)
        # Check the neighbors of the oxygen (other than the inositol carbon)
        for nbr in oxy_atom.GetNeighbors():
            if nbr.GetIdx() == ring_atom_idx:
                continue
            # Look for a phosphorus atom:
            if nbr.GetSymbol() == "P":
                # Check if this phosphorus has at least one double-bonded oxygen.
                dblO_count = 0
                for p_nbr in nbr.GetNeighbors():
                    # Skip oxygen that is the same as our oxy_atom because we came from there.
                    if p_nbr.GetSymbol() != "O":
                        continue
                    bond = mol.GetBondBetweenAtoms(nbr.GetIdx(), p_nbr.GetIdx())
                    if bond is not None and bond.GetBondTypeAsDouble() >= 2:
                        dblO_count += 1
                if dblO_count >= 1:
                    phosphate_found = True
                    break
        if phosphate_found:
            break
    if not phosphate_found:
        return False, "No phosphate group found directly attached to an oxygen substituent of the inositol ring"
    
    # All criteria are met.
    return True, "Contains two acyl chains linked via ester bonds, an inositol ring with sufficient oxygen substituents, and a phosphate attached to the inositol ring"

# Example usage (uncomment for testing):
# test_smiles = "[C@@H]1(C(C(C([C@H](C1O)OP(=O)(O)O)O)O)O)OP(OC[C@](COC(CCCCCCCCCCCCCCCCC)=O)([H])OC(CCCCCCCCCCCCCCC)=O)(O)=O"
# print(is_phosphatidylinositol_phosphate(test_smiles))