"""
Classifies: CHEBI:28868 fatty acid anion
"""
#!/usr/bin/env python3
"""
Classifies: Fatty acid anion
Definition: The conjugate base of a fatty acid, arising from deprotonation 
of the carboxylic acid group of the corresponding fatty acid.
A fatty acid anion is expected to have one terminal deprotonated carboxylate group ([C(=O)[O-]])
attached to an aliphatic (acyclic) chain that forms the dominant portion of the structure.
"""

from rdkit import Chem

def is_fatty_acid_anion(smiles: str):
    """
    Determines if a molecule is a fatty acid anion based on its SMILES string.
    
    Revised heuristic criteria:
      - The molecule must be valid.
      - Must contain exactly one deprotonated carboxylate group, identified by the SMARTS "[CX3](=O)[O-]".
      - The carboxylate carbon must be terminal (i.e. having exactly one neighboring carbon).
      - Starting from the α–carbon (the sole carbon neighbor of the carboxylate group), we traverse
        only along carbon–carbon bonds that are acyclic (i.e. bonds not in any ring) and count the length of
        the longest contiguous chain.
      - The longest acyclic chain from the α–carbon must be at least 3 carbons long.
      - To ensure that the long chain is dominant in the non‐cyclic portion of the molecule,
        we compute the ratio: (longest acyclic chain length) / (total acyclic carbon atoms excluding the carboxylate).
        This ratio must be at least 0.70.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a fatty acid anion, False otherwise.
        str: Reason for the classification decision.
    """
    # Parse the SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the pattern for a deprotonated carboxylate group.
    carboxylate_smarts = "[CX3](=O)[O-]"
    carboxylate_pattern = Chem.MolFromSmarts(carboxylate_smarts)
    carboxylate_matches = mol.GetSubstructMatches(carboxylate_pattern)
    
    if len(carboxylate_matches) == 0:
        return False, "No deprotonated carboxylate group ([C(=O)[O-]]) found"
    if len(carboxylate_matches) > 1:
        return False, f"Found {len(carboxylate_matches)} carboxylate groups; expected exactly one"
    
    # Get the carboxylate carbon (first index of the match)
    match = carboxylate_matches[0]
    carboxyl_carbon_idx = match[0]
    carboxyl_carbon = mol.GetAtomWithIdx(carboxyl_carbon_idx)
    
    # Ensure the carboxylate carbon is terminal (exactly one neighboring carbon)
    carbon_neighbors = [nbr for nbr in carboxyl_carbon.GetNeighbors() if nbr.GetAtomicNum() == 6]
    if len(carbon_neighbors) != 1:
        return False, "Carboxylate group is not terminal (expected exactly one carbon neighbor)"
    
    # The α–carbon is the sole carbon neighbor.
    alpha_carbon = carbon_neighbors[0]
    
    # We relax any restrictions on the α–carbon substituents; instead we focus on the acyclic carbon chain.
    # Define a DFS function that only follows acyclic C–C bonds.
    def dfs_longest_chain(atom, visited):
        visited.add(atom.GetIdx())
        max_length = 1  # Count the current atom
        for bond in atom.GetBonds():
            # Only follow bonds that are not in a ring
            if bond.IsInRing():
                continue
            # Get the neighboring atom
            nbr = bond.GetOtherAtom(atom)
            if nbr.GetAtomicNum() != 6:
                continue
            if nbr.GetIdx() in visited:
                continue
            # Continue the DFS on this neighbor
            chain_length = 1 + dfs_longest_chain(nbr, visited.copy())
            if chain_length > max_length:
                max_length = chain_length
        return max_length

    longest_chain = dfs_longest_chain(alpha_carbon, set())
    if longest_chain < 3:
        return False, f"Aliphatic chain appears too short (chain length = {longest_chain})"
    
    # Calculate total number of acyclic carbon atoms in the molecule (excluding the carboxylate carbon).
    acyclic_carbons = [atom for atom in mol.GetAtoms() 
                       if atom.GetAtomicNum() == 6 and not atom.IsInRing()
                       and atom.GetIdx() != carboxyl_carbon_idx]
    if not acyclic_carbons:
        return False, "No acyclic carbon atoms found aside from carboxylate carbon"
    total_acyclic = len(acyclic_carbons)
    
    chain_ratio = longest_chain / total_acyclic
    if chain_ratio < 0.70:
        return False, ("Longest contiguous acyclic carbon chain from the α–carbon is not dominant in the "
                       f"acyclic portion of the structure (chain ratio = {chain_ratio:.2f})")
    
    return True, ("Contains a terminal deprotonated carboxylate group and a dominant acyclic aliphatic chain, "
                  "consistent with a fatty acid anion")


# Example usage (when running as a script)
if __name__ == "__main__":
    # You can test with one of the provided examples, e.g. hexadecanoate:
    test_smiles = "CCCCCCCCCCCCCCCC([O-])=O"
    is_fa, reason = is_fatty_acid_anion(test_smiles)
    print(is_fa, reason)