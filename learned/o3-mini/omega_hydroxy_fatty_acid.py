"""
Classifies: CHEBI:10615 omega-hydroxy fatty acid
"""
#!/usr/bin/env python
"""
Classifies: omega-hydroxy fatty acid
Definition: Any member of the class of naturally-occurring straight-chain fatty acids n carbon atoms long 
with a carboxyl group at position 1 and a hydroxyl at position n (omega).
"""

from rdkit import Chem

def is_omega_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is an omega-hydroxy fatty acid based on its SMILES string.
    Criteria:
      - The molecule must have a carboxylic acid group (-C(=O)O).
      - The molecule must have a terminal hydroxyl (-OH) group at the chain end opposite 
        to the carboxyl group.
    
    Approach:
      1. Parse the SMILES string and add explicit hydrogens.
      2. Identify the acid group using a specific SMARTS pattern, and record the acid carbon.
      3. Find all carbon atoms (atomic number 6) in the molecule that are terminal on the carbon chain 
         (i.e. have exactly one neighboring carbon).
      4. Exclude the acid carbon from these candidates. The remaining terminal carbon(s) are omega candidates.
      5. For a candidate terminal carbon, check for an attached hydroxyl group (-OH) 
         (an oxygen attached via a single bond with at least one hydrogen).
    
    Args:
        smiles (str): SMILES string of the molecule.
    Returns:
        bool: True if the molecule qualifies as an omega-hydroxy fatty acid, False otherwise.
        str: Explanation with the classification decision.
    """
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens so that hydroxyl hydrogens are represented
    mol = Chem.AddHs(mol)
    
    # Identify the carboxylic acid group with a SMARTS pattern.
    # Using a specific pattern for a proper acid: [CX3](=O)[OX2H1]
    acid_smarts = "[CX3](=O)[OX2H1]"
    acid_pattern = Chem.MolFromSmarts(acid_smarts)
    acid_matches = mol.GetSubstructMatches(acid_pattern)
    if not acid_matches:
        return False, "Missing carboxylic acid group (-C(=O)O)"
    
    # Record all acid carbon indices (the first atom in the match is the carbon for our pattern)
    acid_carbon_idxs = {match[0] for match in acid_matches}
    
    # Identify terminal carbons in the molecule.
    # In a straight-chain fatty acid, terminal carbons on the chain have exactly one carbon neighbor.
    terminal_carbons = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:  # Only consider carbons
            continue
        # Count neighbors that are carbon atoms
        carbon_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 6]
        if len(carbon_neighbors) == 1:
            terminal_carbons.append(atom.GetIdx())
    
    # Expect two terminal carbons: one belongs to the acid group and the other should be the omega carbon.
    omega_candidates = [idx for idx in terminal_carbons if idx not in acid_carbon_idxs]
    if not omega_candidates:
        return False, "No terminal carbon (other than the carboxyl carbon) found in the chain"
    
    # For each candidate terminal carbon, check if it has an attached hydroxyl (-OH) group.
    for idx in omega_candidates:
        carbon = mol.GetAtomWithIdx(idx)
        for nbr in carbon.GetNeighbors():
            # Look for neighboring oxygen atoms
            if nbr.GetAtomicNum() != 8:
                continue
            bond = mol.GetBondBetweenAtoms(carbon.GetIdx(), nbr.GetIdx())
            # Check for a single bond between the candidate carbon and the oxygen.
            if bond is None or bond.GetBondType() != Chem.BondType.SINGLE:
                continue
            # Check if the oxygen has at least one hydrogen (i.e. is an -OH group)
            if nbr.GetTotalNumHs() >= 1:
                return True, ("Found terminal hydroxyl (-OH) group on the terminal carbon opposite to the "
                              "carboxyl group, consistent with an omega-hydroxy fatty acid.")
    
    return False, "No terminal hydroxyl (-OH) group found on the opposite end of the chain from the carboxyl group"

# Example usage (uncomment for testing):
# if __name__ == "__main__":
#     test_smiles = "OCCCCCCCC(O)=O"  # Example: 8-hydroxyoctanoic acid
#     result, reason = is_omega_hydroxy_fatty_acid(test_smiles)
#     print(result, reason)