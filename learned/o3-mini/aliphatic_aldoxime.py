"""
Classifies: CHEBI:82744 aliphatic aldoxime
"""
"""
Classifies: Any aldoxime derived from an aliphatic aldehyde.
Aldoxime = group of the form CH=N–OH where the carbon is derived from an aldehyde.
An aliphatic aldehyde is defined as an aldehyde that is not aromatic and is not substituted by (or connected to) any aromatic group.
The strategy:
  1. Use a SMARTS pattern "[CH1;!a]=[N][O;H]" to find aldoxime moieties,
     which ensures the carbon has one hydrogen and is not flagged as aromatic.
  2. For each match, check:
      – The aldehyde carbon is not aromatic.
      – It has exactly one heavy-atom substituent aside from the nitrogen.
      – The substituent is not directly joined to an aromatic bond or atom.
      – Additionally, check the immediate neighbors of that substituent for any aromatic atoms.
If any match passes all these filters then the molecule will be classified as an aliphatic aldoxime.
"""

from rdkit import Chem

def is_aliphatic_aldoxime(smiles: str):
    """
    Determines if a molecule is an aliphatic aldoxime based on its SMILES string.
    An aliphatic aldoxime is any aldoxime (CH=N–OH) that comes from an aliphatic aldehyde.
    The carbon involved should be a CH unit (with one hydrogen) that is not aromatic and is not 
    substituted by atoms directly connected (by aromatic bonds or belonging to aromatic systems) 
    to aromatic atoms.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an aliphatic aldoxime, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens to correctly check the CH count.
    mol = Chem.AddHs(mol)
    
    # Define a SMARTS for an aldoxime group.
    # The pattern: [CH1;!a]=[N][O;H]
    aldoxime_smarts = "[CH1;!a]=[N][O;H]"
    aldoxime_pattern = Chem.MolFromSmarts(aldoxime_smarts)
    if aldoxime_pattern is None:
        return False, "Error in defining SMARTS for aldoxime"
    
    # Find all matches for the aldoxime group.
    matches = mol.GetSubstructMatches(aldoxime_pattern)
    if not matches:
        return False, "No aldoxime group (CH=N–OH with a CH unit) found"
    
    # Helper function: check if a bond is aromatic.
    def bond_is_aromatic(mol, idx1, idx2):
        bond = mol.GetBondBetweenAtoms(idx1, idx2)
        if bond is not None and bond.GetIsAromatic():
            return True
        return False

    # For each aldoxime match, apply further aliphatic filters.
    for match in matches:
        # match returns indices for (carbon, nitrogen, oxygen) per our SMARTS order.
        carbon_idx, nitrogen_idx, oxygen_idx = match
        carbon_atom = mol.GetAtomWithIdx(carbon_idx)
        
        # 1. The carbon itself must not be aromatic.
        if carbon_atom.GetIsAromatic():
            continue  # not aliphatic
        
        # 2. Check that the aldehyde carbon is only connected to hydrogen(s) and one heavy substituent apart from the nitrogen.
        heavy_neighbors = [n for n in carbon_atom.GetNeighbors() if n.GetAtomicNum() != 1 and n.GetIdx() != nitrogen_idx]
        if len(heavy_neighbors) != 1:
            # Either too many substituents or none.
            continue
        
        # 3. Examine that one heavy neighbor (the rest group) is not directly attached via an aromatic bond.
        neighbor = heavy_neighbors[0]
        if bond_is_aromatic(mol, carbon_idx, neighbor.GetIdx()):
            continue
        if neighbor.GetIsAromatic():
            continue
        
        # 4. Check one level deeper: none of the neighbor's non-hydrogen neighbors (excluding the carbon) should be aromatic
        qualifies = True
        for sub_neighbor in neighbor.GetNeighbors():
            if sub_neighbor.GetIdx() == carbon_idx:
                continue
            if sub_neighbor.GetAtomicNum() == 1:
                continue
            # Also check the bond between neighbor and sub_neighbor
            if bond_is_aromatic(mol, neighbor.GetIdx(), sub_neighbor.GetIdx()):
                qualifies = False
                break
            if sub_neighbor.GetIsAromatic():
                qualifies = False
                break
        if not qualifies:
            continue

        # If we pass all tests, then we consider the aldoxime as derived from an aliphatic aldehyde.
        return True, "Aldoxime group derived from an aliphatic aldehyde found"
    
    return False, "Aldoxime group(s) found but substituents suggest aromatic or conjugated system"

# Optional: Testing examples if run as a script.
if __name__ == "__main__":
    # (Z)-acetaldehyde oxime
    test_smiles = "[H]\\C(C)=N\\O"
    result, reason = is_aliphatic_aldoxime(test_smiles)
    print(result, "->", reason)