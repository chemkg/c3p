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
    Determines if a molecule is a very long-chain fatty acid (VLCFA) based on its SMILES string.
    
    Our approach:
      1. Parse the SMILES string into an RDKit molecule.
      2. Check for the presence of a carboxylic acid group using two SMARTS patterns 
         (protonated and deprotonated).
      3. Verify that the molecule’s heavy-atom composition is predominantly carbon.
      4. For each acid group found, start from the acid carbon and compute the longest acyclic 
         chain length by following only carbon atoms that are not part of any ring.
      5. If the maximum chain length (number of connected carbons, including the acid carbon)
         is greater than 22, the molecule qualifies as a very long-chain fatty acid.
         If the chain length is also greater than 27, we note that it is an ultra-long-chain species.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is a very long-chain fatty acid, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS patterns for the carboxyl group.
    acid_patterns = [
        Chem.MolFromSmarts("[CX3](=O)[OX2H]"),   # protonated carboxylic acid
        Chem.MolFromSmarts("[CX3](=O)[OX1-]")      # deprotonated carboxylate
    ]
    acid_matches = []
    for pattern in acid_patterns:
        acid_matches.extend(mol.GetSubstructMatches(pattern))
    if not acid_matches:
        return False, "No carboxylic acid group found"
    
    # Check overall carbon content of non-hydrogen atoms.
    heavy_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() != 1]
    total_heavy = len(heavy_atoms)
    total_carbons = sum(1 for atom in heavy_atoms if atom.GetAtomicNum() == 6)
    carbon_fraction = total_carbons / total_heavy if total_heavy > 0 else 0.0
    if carbon_fraction < 0.8:
        return False, f"Carbon fraction too low ({carbon_fraction:.2f}); not typical of a fatty acid"
    
    # Define a function to compute the longest acyclic carbon chain 
    # starting from a given atom, following only carbon atoms that are not in any ring.
    def longest_chain(atom, visited):
        # current atom is counted (should be carbon by selection)
        max_length = 1
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6 and (not nbr.IsInRing()) and nbr.GetIdx() not in visited:
                new_visited = visited.copy()
                new_visited.add(nbr.GetIdx())
                length = 1 + longest_chain(nbr, new_visited)
                if length > max_length:
                    max_length = length
        return max_length

    max_chain_len = 0
    # For each acid group match, start from the acid carbon (first atom in SMARTS match).
    for match in acid_matches:
        acid_c = mol.GetAtomWithIdx(match[0])
        # If acid carbon is in a ring, skip it (fatty acids should have an acyclic chain)
        if acid_c.IsInRing():
            continue
        chain_len = longest_chain(acid_c, {acid_c.GetIdx()})
        if chain_len > max_chain_len:
            max_chain_len = chain_len

    # Decision based on the maximum acyclic carbon chain length.
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