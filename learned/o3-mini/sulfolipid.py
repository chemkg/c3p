"""
Classifies: CHEBI:61384 sulfolipid
"""
"""
Classifies: sulfolipid
Definition: A compound containing a sulfonic acid residue (–OS(=O)(=O)O or its deprotonated variant)
            joined by a carbon–sulfur bond to a lipid (here defined as having at least 10 contiguous aliphatic carbons).
            
The strategy is:
  1. Parse the SMILES.
  2. Look for a sulfonic acid group that is attached to a carbon (using two SMARTS patterns).
  3. For each candidate, perform a depth-first search (DFS) to find the longest contiguous chain of non‐aromatic,
     sp3 (aliphatic) carbons connected by single bonds. The chain must contain at least 10 carbons.
  4. Return True if any candidate passes these criteria.
  
If either the sulfo group is missing or no candidate carbon is directly attached to a long aliphatic chain,
the function returns False with a reason.
"""

from rdkit import Chem

def is_sulfolipid(smiles: str):
    """
    Determines if a molecule is a sulfolipid by verifying that:
      (a) a sulfonic acid group (–OS(=O)(=O)O or its deprotonated variant) is present,
      (b) the sulfonic acid group is attached via a carbon–sulfur bond,
      (c) the connecting carbon is part of a long aliphatic chain (>= 10 contiguous sp3 carbons).
      
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a sulfolipid, otherwise False.
        str: Explanation for the decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # Define SMARTS patterns for sulfonic acid group attached to a carbon.
    # We require a carbon directly connected to a sulfur that is in a sulfonic acid (protonated or deprotonated) state.
    # Two patterns: one for deprotonated and one for protonated forms.
    sulfo_deprot = Chem.MolFromSmarts("[#6]-S(=O)(=O)[O-]")
    sulfo_prot   = Chem.MolFromSmarts("[#6]-S(=O)(=O)O")
    
    # Collect candidate matches: each is a tuple (C_idx, S_idx) where the carbon is the first atom.
    candidate_matches = []
    for pattern in (sulfo_deprot, sulfo_prot):
        matches = mol.GetSubstructMatches(pattern)
        for match in matches:
            # We assume the pattern is defined such that the first atom is the carbon.
            candidate_matches.append(match)
    
    if not candidate_matches:
        return False, "No sulfonic acid group attached to a carbon (C–S bond) found."
    
    # For each candidate match, the first atom (index match[0]) is the carbon attached to the sulfonic group.
    # Now check if that carbon is part of a long aliphatic chain.
    # We define a DFS that traverses only contiguous sp3 (non-aromatic) carbons connected by SINGLE bonds.
    
    def dfs(atom, visited):
        """
        Recursively finds the length of the longest path (chain) in the subgraph of allowed carbon atoms.
        Allowed atoms: atomic number 6, NOT aromatic.
        Allowed bonds: SINGLE bonds.
        
        Args:
            atom: current RDKit atom.
            visited: set of atom indices already visited (to avoid cycles).
        Returns:
            int: Maximum chain-length (number of carbon atoms) found starting from this atom.
        """
        max_length = 1  # count self
        for nbr in atom.GetNeighbors():
            bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
            if nbr.GetAtomicNum() == 6 and (not nbr.GetIsAromatic()) \
               and bond.GetBondType() == Chem.BondType.SINGLE \
               and nbr.GetIdx() not in visited:
                new_visited = visited.copy()
                new_visited.add(nbr.GetIdx())
                chain_length = 1 + dfs(nbr, new_visited)
                if chain_length > max_length:
                    max_length = chain_length
        return max_length

    # For each candidate carbon that is directly attached to a sulfo group, measure the length of its carbon chain.
    for match in candidate_matches:
        carbon_idx = match[0]
        carbon_atom = mol.GetAtomWithIdx(carbon_idx)
        # Only consider aliphatic carbons
        if carbon_atom.GetAtomicNum() != 6 or carbon_atom.GetIsAromatic():
            continue
        # Start DFS from this carbon; mark it as visited.
        chain_length = dfs(carbon_atom, {carbon_idx})
        if chain_length >= 10:
            return True, f"Contains sulfonic acid group attached via C–S bond to a lipid chain of {chain_length} carbons."
    
    return False, "No sulfolipid candidate found: sulfonic acid group is not attached to a sufficiently long aliphatic chain."
    
# Optionally, you can test the function with SMILES strings.
# e.g.:
# result, reason = is_sulfolipid("CCCCCCCCCCCCCCC(=O)O[C@H]1...")  # (example SMILES)
# print(result, reason)