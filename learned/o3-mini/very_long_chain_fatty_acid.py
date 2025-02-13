"""
Classifies: CHEBI:27283 very long-chain fatty acid
"""
#!/usr/bin/env python
"""
Classifies: Very long-chain fatty acid (VLCFA)
Definition: A fatty acid which has a chain length greater than C22.
            Very long-chain fatty acids which have a chain length greater than C27 
            are also known as ultra-long-chain fatty acids.
The function is_very_long_chain_fatty_acid takes a SMILES string as input and 
returns a tuple: (bool, reason).
The approach:
  1. Parse the SMILES string.
  2. Look for a carboxylic acid group (either protonated or deprotonated); if missing, fail.
  3. Ensure the overall heavy atoms are predominantly carbon.
  4. For each acid group found, pick the acid (carbonyl) carbon and then select 
     one or more neighboring carbon(s) as the start of the tail.
  5. Use a depth-first search that (a) only follows carbon atoms not in rings and 
     (b) returns the longest connected chain (as a path of atom indices). 
  6. Include the acid carbon in the chain to yield a “chain length.”
  7. Also, compare the length of the longest chain to the total number of carbons in 
     the molecule. For a “pure” fatty acid this ratio should be high (≥0.70).
  8. Based on the chain length we decide if the candidate is a very long-chain fatty acid,
     and we mention if it qualifies as an ultra-long-chain fatty acid.
If any of these checks fails the program returns False plus the reason.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_very_long_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a very long-chain fatty acid based on its SMILES string.
    
    Approach:
      - The molecule must have a carboxylic acid group (protonated or deprotonated).
      - From the acid carbon, follow the connected chain of carbon atoms (ignoring rings)
        to compute the longest chain.
      - In a typical fatty acid, nearly all the carbons are part of this chain. So we require
        that the longest chain accounts for at least 70% of all carbon atoms.
      - Finally, the chain length (counting the acid carbon) must be greater than 22.
        If it exceeds 27, we note that it is an ultra-long-chain fatty acid.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        (bool, str): True/False and a text reason for the classification.
    """
    # Parse the molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS patterns for the carboxyl group:
    # [CX3](=O)[OX2H] covers protonated acid; [CX3](=O)[OX1-] covers deprotonated acid.
    acid_patterns = [
        Chem.MolFromSmarts("[CX3](=O)[OX2H]"),
        Chem.MolFromSmarts("[CX3](=O)[OX1-]")
    ]
    acid_matches = []
    for patt in acid_patterns:
        acid_matches.extend(mol.GetSubstructMatches(patt))
    if not acid_matches:
        return False, "No carboxylic acid group found"
    
    # Ensure the non-hydrogen heavy atoms are mostly carbon:
    heavy_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() != 1]
    total_heavy = len(heavy_atoms)
    total_carbons = sum(1 for atom in heavy_atoms if atom.GetAtomicNum() == 6)
    if total_heavy == 0 or (total_carbons / total_heavy) < 0.8:
        return False, f"Carbon fraction too low ({total_carbons/total_heavy if total_heavy else 0:.2f}); not typical of a fatty acid"

    # Define a DFS function to traverse carbon atoms (not in rings) and find the longest path.
    # It returns a tuple (max_length, path) where path is a list of atom indices.
    def dfs(atom, visited):
        best_len = 1
        best_path = [atom.GetIdx()]
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() != 6:
                continue
            # Only follow atoms not yet visited and (preferably) not in a ring.
            if nbr.GetIdx() in visited or nbr.IsInRing():
                continue
            new_visited = visited | {nbr.GetIdx()}
            l, path = dfs(nbr, new_visited)
            if 1 + l > best_len:
                best_len = 1 + l
                best_path = [atom.GetIdx()] + path
        return best_len, best_path

    candidate_chain_length = 0
    candidate_chain_path = []

    # For each acid match, get the acid carbon (first atom in the SMARTS match).
    # In typical fatty acids the acid carbon has two neighbors: one oxygen and one carbon.
    for match in acid_matches:
        acid_atom = mol.GetAtomWithIdx(match[0])
        # If the acid carbon is in a ring, skip it.
        if acid_atom.IsInRing():
            continue
        # Look for carbon neighbors (tail start)
        for nbr in acid_atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6:
                # Start DFS from this carbon neighbor with visited including acid_atom and neighbor.
                cur_len, cur_path = dfs(nbr, visited={acid_atom.GetIdx(), nbr.GetIdx()})
                # Include the acid carbon in the chain length.
                total_len = cur_len + 1
                full_path = [acid_atom.GetIdx()] + cur_path
                if total_len > candidate_chain_length:
                    candidate_chain_length = total_len
                    candidate_chain_path = full_path

    # If we never found any carbon neighbor to the acid carbon, this is not a fatty acid.
    if candidate_chain_length == 0:
        return False, "Could not determine a carbon chain starting from the acid group"

    # Filter on overall structure:
    # For a typical fatty acid, the main chain should account for most carbon atoms.
    # We require that the ratio of chain carbons to total carbons be at least 0.70.
    chain_ratio = candidate_chain_length / total_carbons
    if chain_ratio < 0.70:
        return False, f"Longest chain accounts for only {chain_ratio:.2f} fraction of carbons; structure appears too complex to be a simple fatty acid"

    # Check chain length with respect to the definition (>22 carbons including the acid carbon).
    if candidate_chain_length <= 22:
        return False, f"Carbon count is {candidate_chain_length}, which is not greater than 22 required for a very long-chain fatty acid"
    
    # Classify as ultra-long if chain length is greater than 27.
    if candidate_chain_length > 27:
        reason = (f"Contains a carboxylic acid group and has {candidate_chain_length} connected carbons "
                  f"(ultra-long-chain fatty acid)")
    else:
        reason = (f"Contains a carboxylic acid group and has {candidate_chain_length} connected carbons, "
                  f"which qualifies as a very long-chain fatty acid")
    
    return True, reason

# Example usage:
if __name__ == "__main__":
    example_smiles = "OC(=O)CCCCC#CCC#CCC#CCC#CCC#CCCCCC"  # 6,9,12,15,18-Tetracosapentaynoic acid
    valid, msg = is_very_long_chain_fatty_acid(example_smiles)
    print(valid, msg)