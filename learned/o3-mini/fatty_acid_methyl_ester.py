"""
Classifies: CHEBI:4986 fatty acid methyl ester
"""
#!/usr/bin/env python3
"""
Classifies: fatty acid methyl ester (FAME)

A fatty acid methyl ester is defined as a fatty acid ester obtained by the formal condensation
of a fatty acid with methanol. The approach used here is:
  1. Look for an ester group matching the SMARTS pattern: [C:1](=[O:2])[O:3][C:4].
  2. Check that the alkoxy fragment ([C:4]) is a methyl group.
  3. From the carbonyl carbon (atom 1), identify the unique acyl chain (the carbon neighbor
     other than the ester oxygen) by using a DFS that does not cross the bond to the ester oxygen.
  4. Verify that the acyl fragment has at least 4 carbon atoms, is acyclic and does not contain
     an extra carboxyl (free acid) group.
Any failure to meet these conditions means the molecule is not recognised as a FAME.
"""
from rdkit import Chem

def is_fatty_acid_methyl_ester(smiles: str):
    """
    Determines if a molecule is a fatty acid methyl ester (FAME) based on its SMILES.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if defined as a FAME, False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Define ester SMARTS: [C:1](=[O:2])[O:3][C:4]
    ester_smarts = "[C:1](=[O:2])[O:3][C:4]"
    ester_pattern = Chem.MolFromSmarts(ester_smarts)
    matches = mol.GetSubstructMatches(ester_pattern)
    if not matches:
        return False, "No ester functional group matching [C](=[O])[O]C found"
    
    # Helper function: check if an atom is a methyl group
    def is_methyl(atom):
        if atom.GetAtomicNum() != 6:
            return False
        # Only one heavy neighbor (typically the ester oxygen)
        heavy_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() != 1]
        return len(heavy_neighbors) == 1

    # Helper function: DFS to collect acyl fragment atom indices
    # Avoid crossing the forbidden bond (the bond connecting carbonyl to ester oxygen)
    def dfs_acyl(atom, forbidden_bond_idx, visited):
        visited.add(atom.GetIdx())
        atoms_set = {atom.GetIdx()}
        for bond in atom.GetBonds():
            # Skip the bond we are not allowed to cross
            if bond.GetIdx() == forbidden_bond_idx:
                continue
            nbr = bond.GetOtherAtom(atom)
            if nbr.GetIdx() in visited:
                continue
            atoms_set.update(dfs_acyl(nbr, forbidden_bond_idx, visited))
        return atoms_set

    # Define a SMARTS pattern for a carboxyl group (free acid), i.e. -C(=O)[OH]
    # We want to ensure no extra carboxyl is present in the acyl fragment.
    free_carboxyl_smarts = Chem.MolFromSmarts("C(=O)[OH]")
    
    # Process each ester match candidate
    for match in matches:
        # In each match:
        # match[0] -> carbonyl carbon, match[1] -> double-bonded oxygen,
        # match[2] -> ester oxygen, match[3] -> alkoxy carbon.
        carbonyl = mol.GetAtomWithIdx(match[0])
        dbl_o = mol.GetAtomWithIdx(match[1])
        ester_oxygen = mol.GetAtomWithIdx(match[2])
        alkoxy = mol.GetAtomWithIdx(match[3])
        
        # Ensure the alkoxy group is a methyl group.
        if not is_methyl(alkoxy):
            continue  # Not a methylester; proceed to next candidate
        
        # Identify the acyl chain: It should be the unique carbon neighbor of the carbonyl (other than the ester oxygen).
        acyl_neighbors = [nbr for nbr in carbonyl.GetNeighbors() 
                          if nbr.GetIdx() != ester_oxygen.GetIdx() and nbr.GetAtomicNum() == 6]
        if len(acyl_neighbors) != 1:
            continue  # Ambiguous or no acyl chain found.
        acyl_start = acyl_neighbors[0]
        
        # Reject if the carbonyl or acyl start are in a ring.
        if carbonyl.IsInRing() or acyl_start.IsInRing():
            continue
        
        # Identify the bond between carbonyl and ester oxygen so that we don't cross it during DFS.
        bond_to_block = mol.GetBondBetweenAtoms(carbonyl.GetIdx(), ester_oxygen.GetIdx())
        if bond_to_block is None:
            continue
        forbidden_bond_idx = bond_to_block.GetIdx()
        
        # Perform DFS starting from the carbonyl to get the entire acyl fragment.
        acyl_atom_indices = dfs_acyl(carbonyl, forbidden_bond_idx, set())
        if not acyl_atom_indices:
            continue  # Failed to retrieve acyl fragment.
        
        # Count carbons in the acyl fragment (should be >= 4, counting the carbonyl carbon).
        n_carbons = sum(1 for idx in acyl_atom_indices if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
        if n_carbons < 4:
            continue  # Acyl fragment too short.
        
        # Ensure none of the atoms in the acyl fragment are in a ring
        if any(mol.GetAtomWithIdx(idx).IsInRing() for idx in acyl_atom_indices):
            continue

        # Check for additional free carboxyl groups within the acyl fragment.
        # We search for the free acid SMARTS in the entire molecule and count only those matches fully within the acyl fragment.
        extra_carboxyl_count = 0
        for carboxyl_match in mol.GetSubstructMatches(free_carboxyl_smarts):
            if set(carboxyl_match).issubset(acyl_atom_indices):
                extra_carboxyl_count += 1
        if extra_carboxyl_count > 0:
            continue  # Presence of a free carboxyl group suggests a diacid ester or a different scaffold.
        
        # If all the conditions are met, classify as FAME.
        return True, "Contains a methoxy ester group with an appropriate acyl chain (FAME identified)"
    
    return False, "Ester group found but no candidate met all FAME criteria"


# (Optional) Quick testing when module is run directly.
if __name__ == "__main__":
    # Test with an example (methyl octanoate)
    test_smiles = "O=C(OC)CCCCCCCC"
    result, reason = is_fatty_acid_methyl_ester(test_smiles)
    print(result, reason)