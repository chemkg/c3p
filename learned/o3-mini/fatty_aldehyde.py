"""
Classifies: CHEBI:35746 fatty aldehyde
"""
"""
Classifies: Fatty Aldehydes
Definition: An aldehyde formally arising from reduction of the carboxylic acid group of its corresponding fatty acid, 
having a carbonyl group at one end of the carbon chain.
This improved classifier uses multiple heuristics:
  - The molecule must be parsed and be acyclic.
  - It must have a minimum number of carbons and a low oxygen content relative to carbons.
  - It must contain at least one aldehyde group ([CX3H1](=O)) where the carbonyl carbon is attached to exactly one carbon.
  - From that carbon neighbor (“alpha carbon”), we compute the longest contiguous carbon chain (even if some branching occurs)
    and require that it meet a minimum length.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_fatty_aldehyde(smiles: str):
    """
    Determines if a molecule is a fatty aldehyde based on its SMILES string.
    A fatty aldehyde should originate from a fatty acid reduction and has a terminal aldehyde group
    attached to a long aliphatic (or nearly linear) chain.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule meets all criteria, False otherwise.
        str: Explanation for the decision.
    """
    
    # Parse molecule and ensure it is valid.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Reject molecules with rings (fatty aldehydes are typically acyclic).
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule contains rings, which is not typical for a fatty aldehyde"
    
    # Heuristic: require a minimum number of carbon atoms.
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 6:
        return False, f"Not enough carbon atoms ({c_count}) for a fatty aldehyde"
    
    # Heuristic: require a minimum molecular weight.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 100:
        return False, f"Molecular weight too low ({mol_wt:.2f} Da) for a fatty aldehyde"
    
    # Additional check: fatty compounds are relatively nonpolar.
    # We use the ratio of carbons to oxygens to help filter out highly oxidized molecules (like sugar acids).
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    # If there are oxygens, require at least 2x more carbons than oxygens.
    if o_count > 0 and (c_count / o_count) < 2:
        return False, f"Carbon-to-oxygen ratio too low (C:{c_count} O:{o_count}) for a fatty aldehyde"
    
    # Look for aldehyde groups with SMARTS pattern.
    aldehyde_pattern = Chem.MolFromSmarts("[CX3H1](=O)")
    matches = mol.GetSubstructMatches(aldehyde_pattern)
    if not matches:
        return False, "No aldehyde group found"
    
    # Define minimum required carboxylic chain length beyond the aldehyde.
    # This is the number of carbon atoms in the longest path starting from the carbon attached to the aldehyde.
    MIN_CHAIN_CARBONS = 4
    
    # Build a simple helper to compute longest simple path (in terms of number of carbons) 
    # in the carbon-only subgraph. The molecule is acyclic so this DFS is safe.
    def dfs_longest_chain(curr_idx, visited):
        curr_atom = mol.GetAtomWithIdx(curr_idx)
        max_length = 1  # count current atom
        for nbr in curr_atom.GetNeighbors():
            # Consider only carbon neighbors that have not been visited.
            if nbr.GetAtomicNum() == 6 and nbr.GetIdx() not in visited:
                length = 1 + dfs_longest_chain(nbr.GetIdx(), visited | {nbr.GetIdx()})
                if length > max_length:
                    max_length = length
        return max_length

    # Now check each aldehyde candidate to see if it is terminal.
    terminal_aldehyde_found = False
    for match in matches:
        # In the SMARTS [CX3H1](=O), the match returns the carbonyl carbon first.
        aldehyde_c = mol.GetAtomWithIdx(match[0])
        # Identify the carbon neighbor(s) (ignore oxygen neighbors)
        neighbor_carbons = [nbr for nbr in aldehyde_c.GetNeighbors() if nbr.GetAtomicNum() == 6]
        # For a terminal aldehyde group, the carbonyl carbon should have exactly one carbon neighbor.
        if len(neighbor_carbons) != 1:
            continue
        alpha = neighbor_carbons[0]
        # Compute the longest carbon chain starting from the alpha carbon.
        chain_length = dfs_longest_chain(alpha.GetIdx(), {alpha.GetIdx()})
        # We require that the chain (starting from the alpha carbon) be long enough.
        if chain_length >= MIN_CHAIN_CARBONS:
            terminal_aldehyde_found = True
            break

    if not terminal_aldehyde_found:
        return (False, "Aldehyde group is present but does not appear terminal or is attached to an insufficiently long aliphatic chain")
    
    return True, "Contains a terminal aldehyde group attached to a long aliphatic chain and qualifies as a fatty aldehyde"

# Example usage (uncomment to test):
# test_smiles = "O=CCCCCCCCCC/C=C/CC"  # 11E-Tetradecenal
# print(is_fatty_aldehyde(test_smiles))