"""
Classifies: CHEBI:27283 very long-chain fatty acid
"""
#!/usr/bin/env python
"""
Classifies: A very long-chain fatty acid
Definition: A fatty acid which has a chain length greater than C22.
            Very long-chain fatty acids which have a chain length greater than C27 are also known as ultra-long-chain fatty acids.
The function is_very_long_chain_fatty_acid takes a SMILES string as input and returns a tuple: (bool, reason)
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_very_long_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a very long-chain fatty acid based on its SMILES string.
    
    Our approach:
      1. Parse the SMILES string into an RDKit molecule.
      2. Check for the presence of a carboxylic acid group.
          We use the SMARTS pattern 'C(=O)[O;H]' – note that in RDKit hydrogens may be implicit,
          so this pattern is a heuristic to ensure a carboxyl function is present.
      3. Count the number of carbon atoms (atomic number 6) as a proxy for the fatty acid chain length.
         (This simple approach assumes that the majority of carbons are in the acyl chain.)
      4. If the total carbon count is greater than 22 then we classify the compound as a very 
         long-chain fatty acid. If the count exceeds 27 we note that this is in the ultra-long-chain range.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a very long-chain fatty acid, False otherwise.
        str: Reason for classification.
    """
    # Parse the SMILES string into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for presence of carboxylic acid group using a SMARTS pattern.
    # This pattern looks for a carbon with a double-bond to oxygen and a single-bond to an OH.
    acid_pattern = Chem.MolFromSmarts("C(=O)[O;H]")
    if not mol.HasSubstructMatch(acid_pattern):
        return False, "No carboxylic acid group found"
    
    # Count total number of carbon atoms in the molecule.
    # (This is a simplification; in a real analysis one might try to extract the longest continuous alkyl chain attached to the COOH.)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    
    # Determine whether the number of carbons is sufficient.
    if c_count <= 22:
        return False, f"Carbon count is {c_count}, which is not greater than 22 required for a very long-chain fatty acid"
    
    # Prepare a message based on whether the molecule qualifies as ultra-long-chain.
    if c_count > 27:
        reason = f"Contains a carboxylic acid group and has {c_count} carbons (ultra-long-chain fatty acid)"
    else:
        reason = f"Contains a carboxylic acid group and has {c_count} carbons, which qualifies as a very long-chain fatty acid"
    
    return True, reason

# Example usage (you can test the function with one of the provided examples):
if __name__ == "__main__":
    test_smiles = "OC(=O)CCCCC#CCC#CCC#CCC#CCC#CCCCCC"  # 6,9,12,15,18-Tetracosapentaynoic acid
    result, message = is_very_long_chain_fatty_acid(test_smiles)
    print(result, message)