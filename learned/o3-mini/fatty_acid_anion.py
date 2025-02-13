"""
Classifies: CHEBI:28868 fatty acid anion
"""
#!/usr/bin/env python3
"""
Classifies: Fatty acid anion
Definition: The conjugate base of a fatty acid, arising from deprotonation 
of the carboxylic acid group of the corresponding fatty acid.
A fatty acid anion is expected to have exactly one terminal deprotonated carboxylate group ([C(=O)[O-]])
attached to an acyclic aliphatic chain that forms the dominant portion of the molecule.
"""

from rdkit import Chem

def is_fatty_acid_anion(smiles: str):
    """
    Determines if a molecule is a fatty acid anion based on its SMILES string.
    
    Revised heuristic criteria:
      - The molecule must be valid.
      - It must contain exactly one deprotonated carboxylate group,
        identified by the SMARTS "[CX3](=O)[O-]".
      - That carboxylate carbon must be terminal, i.e. have exactly one carbon neighbor.
      - From that neighbor (the α–carbon), a depth-first search is performed
        along acyclic (non-ring) C–C bonds to find the longest contiguous carbon chain.
      - The longest chain must be at least 3 carbons long.
      - To ensure that the acyclic chain is “dominant” in the molecule’s non‐cyclic part,
        we compute the ratio: (longest acyclic chain length) / (total non‐cyclic carbons excluding the carboxylate carbon).
        This ratio must be at least 0.50.
      
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a fatty acid anion, False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for the deprotonated carboxylate using the SMARTS "[CX3](=O)[O-]".
    carboxylate_smarts = "[CX3](=O)[O-]"
    carboxylate_pattern = Chem.MolFromSmarts(carboxylate_smarts)
    carboxylate_matches = mol.GetSubstructMatches(carboxylate_pattern)
    
    if len(carboxylate_matches) == 0:
        return False, "No deprotonated carboxylate group ([C(=O)[O-]]) found"
    if len(carboxylate_matches) > 1:
        return False, f"Found {len(carboxylate_matches)} carboxylate groups; expected exactly one"
    
    # The carboxylate match returns a tuple of atom indices; the first index is the carboxyl carbon.
    match = carboxylate_matches[0]
    carboxyl_carbon_idx = match[0]
    carboxyl_carbon = mol.GetAtomWithIdx(carboxyl_carbon_idx)
    
    # Ensure that the carboxylate carbon is terminal by having exactly one neighboring carbon.
    carbon_neighbors = [nbr for nbr in carboxyl_carbon.GetNeighbors() if nbr.GetAtomicNum() == 6]
    if len(carbon_neighbors) != 1:
        return False, "Carboxylate group is not terminal (expected exactly one carbon neighbor)"
    
    # The α–carbon is the sole carbon neighbor.
    alpha_carbon = carbon_neighbors[0]
    
    # Define a recursive depth-first search that explores only along acyclic C–C bonds.
    def dfs_longest_chain(atom, visited):
        visited.add(atom.GetIdx())
        max_length = 1  # count the current atom
        for bond in atom.GetBonds():
            if bond.IsInRing():
                continue  # only follow acyclic bonds
            nbr = bond.GetOtherAtom(atom)
            if nbr.GetAtomicNum() != 6:
                continue  # only follow carbon neighbors
            if nbr.GetIdx() in visited:
                continue
            chain_length = 1 + dfs_longest_chain(nbr, visited.copy())
            if chain_length > max_length:
                max_length = chain_length
        return max_length

    longest_chain = dfs_longest_chain(alpha_carbon, set())
    if longest_chain < 3:
        return False, f"Aliphatic chain appears too short (chain length = {longest_chain})"
    
    # Now, count the total number of non‐cyclic carbon atoms excluding the carboxyl carbon.
    acyclic_carbons = [atom for atom in mol.GetAtoms() 
                       if atom.GetAtomicNum() == 6 and not atom.IsInRing() and atom.GetIdx() != carboxyl_carbon_idx]
    if not acyclic_carbons:
        return False, "No acyclic carbon atoms found aside from the carboxylate carbon"
    total_acyclic = len(acyclic_carbons)
    
    chain_ratio = longest_chain / total_acyclic
    if chain_ratio < 0.50:
        return False, (f"Longest contiguous acyclic carbon chain from the α–carbon is not dominant in the "
                       f"acyclic portion of the structure (chain ratio = {chain_ratio:.2f})")
    
    return True, ("Contains a terminal deprotonated carboxylate group and a dominant acyclic aliphatic chain, "
                  "consistent with a fatty acid anion")


# Example usage
if __name__ == "__main__":
    # Try a known example, e.g. hexadecanoate:
    test_smiles = "CCCCCCCCCCCCCCCC([O-])=O"
    result, explanation = is_fatty_acid_anion(test_smiles)
    print(result, explanation)