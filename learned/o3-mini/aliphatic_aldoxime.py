"""
Classifies: CHEBI:82744 aliphatic aldoxime
"""
"""
Classifies: Any aldoxime derived from an aliphatic aldehyde.
An aldoxime group is defined as CH=N–OH where the corresponding carbon is derived from an aliphatic (non-aromatic)
aldehyde. In particular, the heavy substituent (the R group) on the aldehyde carbon must not contain any aromatic atoms.
The strategy:
  1. Use a SMARTS pattern "[CH1;!a]=[N][O;H]" to identify aldoxime groups,
     which ensures the carbon has one hydrogen and is not flagged as aromatic.
  2. For each match, verify:
      a. The aldehyde carbon is not aromatic.
      b. It is attached to only one heavy substituent (apart from the nitrogen of the oxime).
      c. Recursively check that the substituent group (R group) contains no aromatic atoms.
If any of the matches pass these filters, then we classify the molecule as an aliphatic aldoxime.
"""

from rdkit import Chem

def is_aliphatic_aldoxime(smiles: str):
    """
    Determines if a molecule is an aliphatic aldoxime based on its SMILES string.
    An aliphatic aldoxime is any aldoxime (CH=N–OH) where the aldehyde portion (the CH) is
    derived from a non-aromatic (aliphatic) aldehyde; in other words, its R substituent contains
    no aromatic atoms.

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

    # Add explicit hydrogens to check hydrogen counts correctly.
    mol = Chem.AddHs(mol)

    # Define a SMARTS for an aldoxime group.
    # Pattern: [CH1;!a]=[N][O;H]
    aldoxime_smarts = "[CH1;!a]=[N][O;H]"
    aldoxime_pattern = Chem.MolFromSmarts(aldoxime_smarts)
    if aldoxime_pattern is None:
        return False, "Error in defining SMARTS for aldoxime"
    
    # Find all matches of the aldoxime group.
    matches = mol.GetSubstructMatches(aldoxime_pattern)
    if not matches:
        return False, "No aldoxime group (CH=N–OH with a CH unit) found"
    
    # For each aldoxime match, apply further aliphatic filters.
    for match in matches:
        # According to our SMARTS, match order is (carbon, nitrogen, oxygen)
        carbon_idx, nitrogen_idx, oxygen_idx = match
        carbon_atom = mol.GetAtomWithIdx(carbon_idx)
        
        # 1. The aldehyde carbon must not be aromatic.
        if carbon_atom.GetIsAromatic():
            continue

        # 2. The aldehyde carbon should have exactly one heavy neighbor (non-hydrogen)
        # distinct from the nitrogen atom.
        heavy_neighbors = [nbr for nbr in carbon_atom.GetNeighbors() 
                           if nbr.GetAtomicNum() != 1 and nbr.GetIdx() != nitrogen_idx]
        if len(heavy_neighbors) != 1:
            # Either too many substituents (e.g. substituted by aromatic groups) or none.
            continue
            
        substituent_atom = heavy_neighbors[0]
        
        # 3. Check that the bond directly between the aldehyde carbon and its substituent is not aromatic.
        bond = mol.GetBondBetweenAtoms(carbon_idx, substituent_atom.GetIdx())
        if bond is not None and bond.GetIsAromatic():
            continue
        
        # 4. Traverse the entire substituent (the R group) starting from substituent_atom.
        # Make sure that ONCE we leave the aldoxime carbon, no aromatic atoms are found.
        # We do a breadth-first search (BFS) here.
        visited = set()
        queue = [substituent_atom.GetIdx()]
        aromatic_found = False
        while queue:
            curr_idx = queue.pop(0)
            if curr_idx in visited:
                continue
            visited.add(curr_idx)
            curr_atom = mol.GetAtomWithIdx(curr_idx)
            # If we find an aromatic atom, then the aldehyde is not aliphatic.
            if curr_atom.GetIsAromatic():
                aromatic_found = True
                break
            # Add neighbors that are not the original aldehyde carbon.
            for nbr in curr_atom.GetNeighbors():
                if nbr.GetIdx() == carbon_idx:  # avoid going back to the aldoxime carbon
                    continue
                if nbr.GetAtomicNum() == 1:  # skip hydrogens
                    continue
                if nbr.GetIdx() not in visited:
                    queue.append(nbr.GetIdx())
        if aromatic_found:
            continue

        # Passed all tests: the aldoxime group is attached to an aliphatic residue.
        return True, "Aldoxime group derived from an aliphatic aldehyde found"
    
    # If no aldoxime group passes the tests.
    return False, "Aldoxime group(s) found but substituents suggest aromatic or conjugated system"

# Optional: testing examples when running as a script.
if __name__ == "__main__":
    # Example test: (Z)-acetaldehyde oxime
    test_smiles = "[H]\\C(C)=N\\O"
    result, reason = is_aliphatic_aldoxime(test_smiles)
    print(result, "->", reason)