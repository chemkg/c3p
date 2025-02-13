"""
Classifies: CHEBI:35819 branched-chain fatty acid
"""
"""
Classifies: branched-chain fatty acid.
A branched-chain fatty acid is defined as a fatty acid (carboxylic acid with a long alkyl chain)
that has one or more alkyl substituents branching off from its main (parent) chain.
This implementation checks for the carboxyl group, extracts the carbon chain starting from the
alpha carbon (the carbon directly attached to the carboxyl carbon), and then determines if any
carbon in the main chain has an extra (non‐chain) carbon attached.
"""

from rdkit import Chem

def is_branched_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a branched-chain fatty acid based on its SMILES string.
    
    The algorithm works as follows:
    1. Check that the molecule contains a carboxyl group (the signature of a fatty acid);
       the SMARTS used here will match either a protonated or deprotonated carboxyl.
    2. From the carboxyl group, identify the alpha carbon (the carbon attached to the carboxyl carbon).
       If no such carbon is found then it is not a fatty acid.
    3. Build a subgraph of the molecule containing only carbon atoms. Starting from the alpha carbon,
       we perform a depth-first search (DFS) to find the longest continuous carbon chain. This is our
       candidate fatty acyl chain.
    4. Inspect each carbon in the longest chain: if any carbon has a neighbor that is a carbon not part
       of the main chain (with the exception of the carboxyl carbon attached at the head), then we assume
       a branching substituent is present.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule is a branched-chain fatty acid, False otherwise.
        str: Reason for classification.
    """
    # Parse the molecule from the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Define a SMARTS pattern for carboxylic acid.
    # This pattern matches either the -COOH or -COO[-] groups.
    carboxyl_smarts = "[CX3](=O)[OX1H0-,OX2H1]"
    carboxyl_query = Chem.MolFromSmarts(carboxyl_smarts)
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_query)
    if not carboxyl_matches:
        return False, "No carboxyl group found; not a fatty acid"
    
    # Assume the first match corresponds to the carboxyl group.
    # The SMARTS is constructed so atom 0 is the carbonyl carbon.
    carboxyl_atoms = carboxyl_matches[0]
    carboxyl_C = carboxyl_atoms[0]
    
    # Find the carbon neighbor (alpha carbon) attached to the carboxyl carbon.
    alpha_carbon = None
    carboxyl_atom = mol.GetAtomWithIdx(carboxyl_C)
    for neighbor in carboxyl_atom.GetNeighbors():
        if neighbor.GetAtomicNum() == 6:  # carbon
            alpha_carbon = neighbor.GetIdx()
            break
    if alpha_carbon is None:
        return False, "No alpha carbon attached to carboxyl group; not a fatty acid"
    
    # Build a graph of carbon atoms only.
    carbon_indices = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    carbon_neighbors = {idx: [] for idx in carbon_indices}
    for idx in carbon_indices:
        atom = mol.GetAtomWithIdx(idx)
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6:
                carbon_neighbors[idx].append(nbr.GetIdx())
                
    # Define DFS to find the longest chain (as a list of atom indices) starting from start_idx.
    # We avoid revisiting atoms to prevent cycles.
    def dfs(current, path):
        longest = path
        for nbr in carbon_neighbors[current]:
            if nbr not in path:
                new_path = dfs(nbr, path + [nbr])
                if len(new_path) > len(longest):
                    longest = new_path
        return longest

    # We define the main chain as the longest chain starting from the alpha carbon.
    main_chain = dfs(alpha_carbon, [alpha_carbon])
    
    # To be considered a fatty acid, the chain should be long enough.
    # We expect a fatty acid to have at least 4 carbons in its chain (this threshold can be adjusted).
    if len(main_chain) < 4:
        return False, f"Alkyl chain too short (length {len(main_chain)}) to be a fatty acid"
    
    # Check for branching along the main chain:
    # For each carbon in the main chain, look at its carbon neighbors.
    # Allow the connection to its chain neighbors and, for the alpha carbon, the connection to the carboxyl carbon.
    # Any extra carbon connection is taken as a branch.
    for idx in main_chain:
        atom = mol.GetAtomWithIdx(idx)
        # Count carbon neighbors that are not part of the main chain (or not the carboxyl carbon for the alpha carbon)
        extra_branches = 0
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6:
                nbr_idx = nbr.GetIdx()
                # If this neighbor is not in the main chain or (for the alpha carbon) is not the carboxyl carbon:
                if idx == main_chain[0] and nbr_idx == carboxyl_C:
                    continue  # allow carboxyl connection for the first carbon in chain
                if nbr_idx not in main_chain:
                    extra_branches += 1
        if extra_branches > 0:
            return True, f"Branching detected at carbon atom index {idx}"
            
    # No branching substituents were found along the main alkyl chain.
    return False, "No branching substituents detected on the fatty acyl chain"

# Example usage (you may uncomment and test with one of the provided SMILES strings):
# result, reason = is_branched_chain_fatty_acid("CCC(C)CCCCCCCCCCCCCCCCCCCCCCCCCCC(O)=O")
# print(result, reason)