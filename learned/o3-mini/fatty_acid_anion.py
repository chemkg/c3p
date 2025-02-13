"""
Classifies: CHEBI:28868 fatty acid anion
"""
#!/usr/bin/env python3
"""
Classifies: Fatty acid anion
Definition: The conjugate base of a fatty acid, arising from deprotonation 
of the carboxylic acid group of the corresponding fatty acid.
A fatty acid anion is expected to have one terminal deprotonated carboxylate group ([C(=O)[O-]])
attached to an aliphatic chain that is relatively long and dominant in the molecule.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_fatty_acid_anion(smiles: str):
    """
    Determines if a molecule is a fatty acid anion based on its SMILES string.
    
    Heuristic criteria:
      - Must be a valid molecule.
      - Exactly one deprotonated carboxylate group, defined by the SMARTS "[CX3](=O)[O-]".
      - The carboxylate carbon must be terminal (attached to exactly one carbon).
      - The α–carbon (the carbon attached to the carboxylate carbon) should not lie in a ring
        and should not be overly substituted with atoms other than C or O.
      - A DFS from the α–carbon (only along carbon–carbon bonds and not going back to the carboxylate)
        should yield a contiguous chain whose length is at least 3.
      - In addition the longest chain should represent a “dominant” part of the structure; here we require
        that the ratio of chain-length/(total carbons – 1) is at least 0.7.
      
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if the molecule is classified as a fatty acid anion, False otherwise.
        str: A reason for the classification.
    """
    # Parse the SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS pattern for a deprotonated carboxylate group
    carboxylate_smarts = "[CX3](=O)[O-]"
    carboxylate_pattern = Chem.MolFromSmarts(carboxylate_smarts)
    carboxylate_matches = mol.GetSubstructMatches(carboxylate_pattern)
    
    if len(carboxylate_matches) == 0:
        return False, "No deprotonated carboxylate group ([C(=O)[O-]]) found"
    if len(carboxylate_matches) > 1:
        return False, f"Found {len(carboxylate_matches)} carboxylate groups; expected exactly one"
    
    # The first atom in the match (index 0) is the carbon of the carboxylate group
    match = carboxylate_matches[0]
    carboxyl_carbon_idx = match[0]
    carboxyl_carbon = mol.GetAtomWithIdx(carboxyl_carbon_idx)
    
    # Enforce that the carboxylate carbon is terminal: (should have exactly one carbon neighbor)
    carbon_neighbors = [nbr for nbr in carboxyl_carbon.GetNeighbors() if nbr.GetAtomicNum() == 6]
    if len(carbon_neighbors) != 1:
        return False, "Carboxylate group is not terminal (expected exactly one carbon neighbor)"
    
    alpha_carbon = carbon_neighbors[0]
    
    # Check that the α–carbon is not part of a ring.
    if alpha_carbon.IsInRing():
        return False, "The α–carbon (attached to the carboxylate group) is in a ring"
    
    # Check the direct neighbors of the α-carbon.
    # Allow oxygen substituents (common for hydroxy fatty acids) but disallow atoms like N, S, P, halogens
    forbidden_elements = {7, 15, 16, 17, 35, 53}  # N, P, S, Cl, Br, I (atomic numbers)
    for nbr in alpha_carbon.GetNeighbors():
        if nbr.GetIdx() == carboxyl_carbon.GetIdx():
            continue
        if nbr.GetAtomicNum() not in (6, 8):  # allow C and O
            return False, "The α–carbon is decorated with atoms not typical of a simple aliphatic chain"
    
    # Count total carbon atoms in the molecule.
    total_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if total_carbons < 5:
        return False, f"Too few carbon atoms ({total_carbons}) to be a fatty acid"
    
    # Define a helper that returns the length of the longest contiguous carbon chain
    # starting at a given carbon. We only traverse bonds between carbons.
    def get_longest_chain(atom, visited):
        visited.add(atom.GetIdx())
        max_length = 1  # count current atom
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6 and nbr.GetIdx() not in visited:
                # Do not traverse back to the carboxylate carbon.
                chain_length = 1 + get_longest_chain(nbr, visited.copy())
                if chain_length > max_length:
                    max_length = chain_length
        return max_length

    # Start DFS from the α–carbon.
    chain_length = get_longest_chain(alpha_carbon, set())
    # We require a minimum chain length of 3.
    if chain_length < 3:
        return False, f"Aliphatic chain appears too short (chain length = {chain_length})"
    
    # Heuristic: the longest chain from the α–carbon should form the bulk of the molecule.
    # We subtract the carboxylate carbon from total carbons (as it is “outside” the chain).
    chain_ratio = chain_length / (total_carbons - 1)
    if chain_ratio < 0.70:
        return False, ("Longest contiguous carbon chain from the α–carbon is not dominant "
                       f"in the molecule (chain ratio = {chain_ratio:.2f})")
    
    return True, "Contains a terminal carboxylate group linked to a sufficiently dominant aliphatic chain, consistent with a fatty acid anion"


# Example usage (when running as a script)
if __name__ == "__main__":
    # Try one of the provided examples, e.g. hexadecanoate
    test_smiles = "CCCCCCCCCCCCCCCC([O-])=O"
    is_fa, reason = is_fatty_acid_anion(test_smiles)
    print(is_fa, reason)