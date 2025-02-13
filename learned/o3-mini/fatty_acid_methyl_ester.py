"""
Classifies: CHEBI:4986 fatty acid methyl ester
"""
#!/usr/bin/env python3
"""
Classifies: fatty acid methyl ester (FAME)

A fatty acid methyl ester is defined here as:
  a fatty acid ester that is the carboxylic ester obtained by the formal condensation 
  of a fatty acid with methanol.

Our revised approach looks for an ester group matching the SMARTS [C:1](=[O:2])[O:3][C:4],
ensures that the alkoxy group attached at atom :4 is a methyl group, and then
focuses on the acyl chain connected to the carbonyl carbon. We require that:
  • the carbonyl carbon has exactly one carbon substituent aside from the ester oxygen,
  • the fatty acid substituent (i.e. acyl chain) is not part of any ring (thus is acyclic),
  • and the acyl chain contains at least three carbon atoms.
These criteria help distinguish FAMEs from other esters that might contain a methoxy group.
"""

from rdkit import Chem

def is_fatty_acid_methyl_ester(smiles: str):
    """
    Determines if a molecule is a fatty acid methyl ester (FAME) based on its SMILES.
    
    The function first searches for an ester substructure defined by the SMARTS
    [C:1](=[O:2])[O:3][C:4] and then:
      - verifies that the alkoxy substituent, atom :4, is a methyl group,
      - checks that the carbonyl carbon (atom :1) has exactly one carbon substituent other than 
        the ester oxygen,
      - and confirms that the acyl chain (the fatty acid residue) is acyclic and contains at least 
        three carbon atoms.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if it is classified as a fatty acid methyl ester, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the SMARTS for an ester group: carbonyl carbon, carbonyl O, ester O, and alkoxy C.
    ester_smarts = "[C:1](=[O:2])[O:3][C:4]"
    ester_pattern = Chem.MolFromSmarts(ester_smarts)
    matches = mol.GetSubstructMatches(ester_pattern)
    if not matches:
        return False, "Ester functional group with methoxy substituent not found"

    # Helper function to check if a given carbon atom is a methyl group:
    # it must be a carbon and only have one heavy (non-hydrogen) neighbor.
    def is_methyl(carbon):
        if carbon.GetAtomicNum() != 6:
            return False
        # Count heavy neighbors
        heavy_neighbors = [nbr for nbr in carbon.GetNeighbors() if nbr.GetAtomicNum() != 1]
        if len(heavy_neighbors) != 1:
            return False
        return True

    # Helper function: perform a DFS over connected carbon atoms (ignoring rings)
    # If any carbon encountered is part of a ring, we mark the chain as invalid.
    def get_acyl_chain_carbons(atom, visited):
        # If this atom is in a ring, return None (invalid chain)
        if atom.IsInRing():
            return None
        # Only count carbon atoms
        if atom.GetAtomicNum() != 6:
            return set()
        chain_carbons = {atom.GetIdx()}
        visited.add(atom.GetIdx())
        for nbr in atom.GetNeighbors():
            # Only continue through carbon atoms that have not been visited.
            if nbr.GetAtomicNum() == 6 and nbr.GetIdx() not in visited:
                # If the neighbor is in a ring, the chain is not a free acyclic fatty acid chain.
                if nbr.IsInRing():
                    return None
                result = get_acyl_chain_carbons(nbr, visited)
                if result is None:
                    return None
                chain_carbons.update(result)
        return chain_carbons

    # Loop over each ester group found.
    for match in matches:
        # match indices: 
        # match[0] -> carbonyl carbon, match[1] -> carbonyl oxygen,
        # match[2] -> ester oxygen, match[3] -> alkoxy carbon.
        carbonyl = mol.GetAtomWithIdx(match[0])
        ester_oxygen = mol.GetAtomWithIdx(match[2])
        alkoxy = mol.GetAtomWithIdx(match[3])

        # Check that the alkoxy group is a methyl group.
        if not is_methyl(alkoxy):
            continue  # not a methoxy group, try next match

        # Identify the acyl chain.
        # Get neighbors of carbonyl that are carbons AND are not the ester oxygen.
        acyl_neighbors = []
        for nbr in carbonyl.GetNeighbors():
            if nbr.GetIdx() == ester_oxygen.GetIdx():
                continue
            if nbr.GetAtomicNum() == 6:
                acyl_neighbors.append(nbr)
        # For a typical fatty acid, there should be exactly one acyl chain attached.
        if len(acyl_neighbors) != 1:
            continue  # ambiguous or missing acyl chain

        acyl_start = acyl_neighbors[0]
        # Reject if the acyl start itself is in a ring.
        if acyl_start.IsInRing():
            continue

        # Use DFS to gather all connected carbon atoms in the acyl chain.
        acyl_chain = get_acyl_chain_carbons(acyl_start, set())
        if acyl_chain is None:
            continue  # part of the chain is cyclic; not a fatty acid chain
        # Require a minimum chain length of 3 carbon atoms (excluding the carbonyl carbon)
        if len(acyl_chain) < 3:
            continue  # acyl chain too short
            
        # If we reach here, the match fulfills our criteria.
        return True, "Contains a methoxy ester group with an acyclic acyl chain of sufficient length (FAME identified)"
    
    # If no matching ester fulfills the criteria, then return False.
    return False, "Ester functional group found but does not match fatty acid methyl ester criteria"

# (Optional) Quick testing when the module is run directly.
if __name__ == "__main__":
    # Example: methyl octanoate should be classified as a FAME.
    test_smiles = "O=C(OC)CCCCCCCC"
    result, reason = is_fatty_acid_methyl_ester(test_smiles)
    print(result, reason)