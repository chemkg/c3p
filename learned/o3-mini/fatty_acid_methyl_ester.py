"""
Classifies: CHEBI:4986 fatty acid methyl ester
"""
#!/usr/bin/env python3
"""
Classifies: fatty acid methyl ester (FAME)

A fatty acid methyl ester is defined here as:
  a fatty acid ester that is the carboxylic ester obtained by the formal condensation 
  of a fatty acid with methanol.
  
This function looks for an ester group of the form [C](=[O])[O][C] in which the [C] 
after the oxygen is a methyl group, and then checks that the acyl part (the fatty acid chain)
attached to the carbonyl carbon is of sufficient length.
"""

from rdkit import Chem

def is_fatty_acid_methyl_ester(smiles: str):
    """
    Determines if a molecule is a fatty acid methyl ester based on its SMILES string.
    
    The function looks for an ester substructure of the form [C:1](=[O:2])[O:3][C:4],
    where the atom :4 (the alkoxy substituent) is a methyl group, and the fatty acid (R group)
    connected to atom :1 (the carbonyl carbon) is of at least minimal alkyl length.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule is a fatty acid methyl ester, False otherwise.
        str: Reason for classification.
    """
    # Parse SMILES to a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS for an ester group: [C:1](=[O:2])[O:3][C:4]
    ester_smarts = "[C:1](=[O:2])[O:3][C:4]"
    ester_pattern = Chem.MolFromSmarts(ester_smarts)
    matches = mol.GetSubstructMatches(ester_pattern)
    if not matches:
        return False, "Ester functional group with methoxy substituent not found"
    
    # Helper function to determine if an atom represents a methyl group 
    # (i.e. a carbon attached only to one heavy atom, so only implicit hydrogens otherwise)
    def is_methyl(carbon):
        # Must be carbon
        if carbon.GetAtomicNum() != 6:
            return False
        # Degree in heavy atoms should be 1 (only attached to the oxygen in the ester)
        if carbon.GetDegree() != 1:
            return False
        return True
    
    # Helper function: DFS to count the number of carbon atoms in the acyl chain
    def count_aliphatic_carbons(atom, visited):
        count = 0
        # Only count if it's a carbon atom (atomic number 6)
        if atom.GetAtomicNum() == 6:
            count = 1
        visited.add(atom.GetIdx())
        for neighbor in atom.GetNeighbors():
            # Continue only if neighbor is carbon, not visited and not the ester oxygen leaving group.
            if neighbor.GetIdx() not in visited and neighbor.GetAtomicNum() == 6:
                count += count_aliphatic_carbons(neighbor, visited)
        return count
    
    # Loop over each ester substructure match; at least one must satisfy our criteria.
    for match in matches:
        # match is a tuple with indices corresponding to the labeled atoms in the SMARTS:
        # match[0]: carbonyl carbon, match[1]: carbonyl oxygen, match[2]: ester oxygen, match[3]: alkoxy carbon.
        carbonyl_carbon = mol.GetAtomWithIdx(match[0])
        ester_oxygen = mol.GetAtomWithIdx(match[2])
        alkoxy_carbon = mol.GetAtomWithIdx(match[3])
        
        # Check that the alkoxy group is a methyl group.
        if not is_methyl(alkoxy_carbon):
            continue  # try next match if available

        # Identify the acyl chain: the fatty acid residue should be attached to the carbonyl carbon
        # aside from the ester oxygen (which is already accounted for). In a carboxylic acid,
        # the carbonyl carbon is bound to exactly three groups: double-bond O, ester oxygen, and the R group.
        acyl_neighbors = [nbr for nbr in carbonyl_carbon.GetNeighbors() 
                          if nbr.GetIdx() not in (ester_oxygen.GetIdx(),) and nbr.GetAtomicNum() == 6]
        if not acyl_neighbors:
            continue  # no fatty acyl chain found in this match
        # Usually there should only be one such neighbor.
        fatty_acid_start = acyl_neighbors[0]
        
        # Count the number of connected carbon atoms in the fatty acid chain.
        # We avoid counting the carbonyl carbon by starting from the neighbor.
        chain_length = count_aliphatic_carbons(fatty_acid_start, set())
        # We require a minimum chain length (e.g., at least 3 carbons in the alkyl part)
        if chain_length < 3:
            continue  # chain too short to be classified as a fatty acid

        # If we reach here, both requirements are satisfied.
        return True, "Contains a methoxy ester group with an acyl chain of sufficient length (FAME identified)"

    # If no match satisfies the criteria:
    return False, "Ester group present but does not match the fatty acid methyl ester criteria"

# (Optional) For quick testing:
if __name__ == "__main__":
    # Example: methyl octanoate: SMILES "O=C(OC)CCCCCCCC" should be True.
    test_smiles = "O=C(OC)CCCCCCCC"
    result, reason = is_fatty_acid_methyl_ester(test_smiles)
    print(result, reason)