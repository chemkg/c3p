"""
Classifies: CHEBI:35746 fatty aldehyde
"""
"""
Classifies: Fatty Aldehydes
Definition: An aldehyde formally arising from reduction of the carboxylic acid group of its corresponding fatty acid, 
having a carbonyl group at one end of the carbon chain.
This classifier applies multiple heuristics:
  - The molecule must be parsed and be acyclic.
  - It must contain a minimum total number of carbons and exceed a minimum molecular weight.
  - The overall heavy-atom composition should be heavily carbon dominated.
  - It must contain at least one terminal aldehyde group ([CX3H1](=O)) where the carbonyl carbon is attached to exactly one carbon.
  - The longest continuous carbon chain (found via DFS over the carbon subgraph) must be “long” (at least 8 carbons).
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_fatty_aldehyde(smiles: str):
    """
    Determines if a molecule is a fatty aldehyde based on its SMILES string.
    A fatty aldehyde is expected to originate from a fatty acid reduction and 
    have a terminal aldehyde group attached to a long, predominantly aliphatic chain.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule meets all criteria, False otherwise.
        str: Explanation for the decision.
    """
    # Parse molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Reject molecules with rings (fatty aldehydes are typically acyclic)
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule contains rings, which is not typical for a fatty aldehyde"

    # Count carbon atoms and total heavy atoms (exclude hydrogens)
    atoms = list(mol.GetAtoms())
    c_count = sum(1 for atom in atoms if atom.GetAtomicNum() == 6)
    heavy_count = sum(1 for atom in atoms if atom.GetAtomicNum() > 1)
    
    # Require a minimum number of carbons (e.g. at least 7)
    if c_count < 7:
        return False, f"Not enough carbon atoms ({c_count}) for a fatty aldehyde"
    
    # Calculate molecular weight and require a minimum threshold (e.g. >= 110 Da)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 110:
        return False, f"Molecular weight too low ({mol_wt:.2f} Da) for a fatty aldehyde"
    
    # Check that the molecule is mostly aliphatic: require a high carbon fraction
    # (For a proper fatty chain, most heavy atoms should be carbons)
    carbon_fraction = c_count / heavy_count if heavy_count else 0
    if carbon_fraction < 0.75:
        return False, f"Carbon fraction too low (C:{c_count} / heavy atoms:{heavy_count} = {carbon_fraction:.2f}) for a fatty aldehyde"
    
    # SMARTS for an aldehyde group where the carbonyl carbon has one bonded hydrogen: [CX3H1](=O)
    aldehyde_pattern = Chem.MolFromSmarts("[CX3H1](=O)")
    aldehyde_matches = mol.GetSubstructMatches(aldehyde_pattern)
    if not aldehyde_matches:
        return False, "No aldehyde group found"
    
    # Define a helper that performs DFS over the carbon-only subgraph to obtain the longest chain length.
    def dfs_longest_chain(curr_idx, visited):
        curr_atom = mol.GetAtomWithIdx(curr_idx)
        max_length = 1  # current atom counts as one
        for nbr in curr_atom.GetNeighbors():
            # Consider only carbon neighbors not yet visited.
            if nbr.GetAtomicNum() == 6 and nbr.GetIdx() not in visited:
                length = 1 + dfs_longest_chain(nbr.GetIdx(), visited | {nbr.GetIdx()})
                if length > max_length:
                    max_length = length
        return max_length

    # Compute overall longest continuous carbon chain in the molecule:
    overall_longest_chain = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:
            chain_length = dfs_longest_chain(atom.GetIdx(), {atom.GetIdx()})
            if chain_length > overall_longest_chain:
                overall_longest_chain = chain_length
    if overall_longest_chain < 8:
        return False, f"Longest continuous carbon chain is too short ({overall_longest_chain} carbons)"

    # Check for at least one terminal aldehyde group.
    # Terminal aldehyde means the carbonyl carbon in the aldehyde group is attached to exactly one carbon neighbor.
    terminal_aldehyde_found = False
    for match in aldehyde_matches:
        aldehyde_c = mol.GetAtomWithIdx(match[0])
        # Get only carbon neighbors (ignore oxygens or hydrogens that are implicit)
        carbon_neighbors = [nbr for nbr in aldehyde_c.GetNeighbors() if nbr.GetAtomicNum() == 6]
        if len(carbon_neighbors) == 1:
            # Optionally, we can check that starting from the alpha carbon (the single neighbor)
            # the longest chain is of sufficient length. Here we require at least 4 carbons extending from it.
            MIN_CHAIN_CARBONS = 4
            alpha = carbon_neighbors[0]
            chain_length = dfs_longest_chain(alpha.GetIdx(), {alpha.GetIdx()})
            if chain_length >= MIN_CHAIN_CARBONS:
                terminal_aldehyde_found = True
                break
    if not terminal_aldehyde_found:
        return False, "Aldehyde group found, but none appears terminal with a sufficiently long alkyl chain"

    return True, "Contains a terminal aldehyde group attached to a sufficiently long, predominantly aliphatic chain and qualifies as a fatty aldehyde"

# Example usage (uncomment the lines below to test):
# test_smiles = "O=CCCCCCCCCC/C=C/CC"  # 11E-Tetradecenal
# result, reason = is_fatty_aldehyde(test_smiles)
# print(result, reason)