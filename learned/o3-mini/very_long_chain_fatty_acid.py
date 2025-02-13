"""
Classifies: CHEBI:27283 very long-chain fatty acid
"""
#!/usr/bin/env python
"""
Classifies: A very long-chain fatty acid
Definition: A fatty acid which has a chain length greater than C22.
            Very long-chain fatty acids which have a chain length greater than C27 
            are also known as ultra-long-chain fatty acids.
The function is_very_long_chain_fatty_acid takes a SMILES string as input and 
returns a tuple: (bool, reason)
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_very_long_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a very long-chain fatty acid based on its SMILES string.
    
    Our approach:
      1. Parse the SMILES string into an RDKit molecule.
      2. Check for the presence of a carboxylic acid group. We use two SMARTS patterns,
         one for the protonated acid ([CX3](=O)[OX2H]) and one for the deprotonated form ([CX3](=O)[OX1-]).
      3. Check that the molecule has a high fraction of carbons among all heavy atoms.
         Fatty acids are composed mostly of carbons (and hydrogens) so we expect a high C : non-H ratio.
      4. For each carboxyl group match, we locate the acid carbon and perform a depth-first search 
         (restricted to carbon atoms only) to extract the connected chain.
         (This is a proxy for the “fatty acyl chain” length; note that in many fatty acids the acid carbon is counted.)
      5. We select the largest such chain and if its length is greater than 22 we classify the molecule
         as a very long-chain fatty acid. If the length is >27 we note it as ultra-long-chain.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a very long-chain fatty acid, False otherwise.
        str: Reason for classification.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS patterns for a carboxyl group (protonated and deprotonated versions).
    acid_patterns = [Chem.MolFromSmarts("[CX3](=O)[OX2H]"),
                     Chem.MolFromSmarts("[CX3](=O)[OX1-]")]

    acid_matches = []
    for pat in acid_patterns:
        acid_matches.extend(mol.GetSubstructMatches(pat))
    if not acid_matches:
        return False, "No carboxylic acid group found"

    # Check overall carbon ratio (non-hydrogen atoms) as a safeguard.
    heavy_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() != 1]
    total_heavy = len(heavy_atoms)
    total_carbons = sum(1 for atom in heavy_atoms if atom.GetAtomicNum() == 6)
    carbon_fraction = total_carbons / total_heavy if total_heavy else 0
    # We expect fatty acids to be mostly carbons; if not, abort.
    if carbon_fraction < 0.8:
        return False, f"Carbon fraction too low ({carbon_fraction:.2f}); not typical of a fatty acid"
    
    # For each acid match, we compute the size of the connected carbon chain from the acid carbon.
    # In the SMARTS the first atom should be the acid (C in C(=O)O); we do a DFS on carbon atoms (atomic num 6).
    def dfs(start_idx, visited):
        visited.add(start_idx)
        count = 1
        for nbr in mol.GetAtomWithIdx(start_idx).GetNeighbors():
            if nbr.GetAtomicNum() == 6 and nbr.GetIdx() not in visited:
                count += dfs(nbr.GetIdx(), visited)
        return count

    max_chain_len = 0
    for match in acid_matches:
        acid_c = match[0]  # the acid carbon atom from the pattern
        # Do DFS starting from acid_c but do not cross over non-carbon atoms.
        visited = set()
        chain_len = dfs(acid_c, visited)
        if chain_len > max_chain_len:
            max_chain_len = chain_len

    # Compare the chain length (number of carbons in the connected acyl chain) with the requirement.
    if max_chain_len <= 22:
        return False, f"Carbon count is {max_chain_len}, which is not greater than 22 required for a very long-chain fatty acid"
    
    if max_chain_len > 27:
        reason = f"Contains a carboxylic acid group and has {max_chain_len} connected carbons (ultra-long-chain fatty acid)"
    else:
        reason = f"Contains a carboxylic acid group and has {max_chain_len} connected carbons, which qualifies as a very long-chain fatty acid"
    
    return True, reason

# Example usage:
if __name__ == "__main__":
    test_smiles = "OC(=O)CCCCC#CCC#CCC#CCC#CCC#CCCCCC"  # Example: 6,9,12,15,18-Tetracosapentaynoic acid
    result, message = is_very_long_chain_fatty_acid(test_smiles)
    print(result, message)