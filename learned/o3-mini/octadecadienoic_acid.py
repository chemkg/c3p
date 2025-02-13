"""
Classifies: CHEBI:25627 octadecadienoic acid
"""
"""
Classifies: Octadecadienoic acid
Definition: Any straight-chain, C18 polyunsaturated fatty acid having exactly two C=C bonds.
The algorithm:
  1. Parses the molecule from SMILES.
  2. Locates all carboxylic acid groups (using two common SMARTS).
  3. For each acid group found:
      • Identify the acid carbon (the matched carbon atom).
      • Find candidate terminal methyl carbons (CH3 groups, i.e. carbons with exactly one carbon neighbor).
      • For each candidate, use a DFS that only follows single/double bonds (excluding ring bonds) and only through carbon atoms.
      • Check that the path length (number of carbons) is exactly 18.
      • For each candidate chain, ensure the chain is “straight”: no extra carbon neighbours (branching) are present.
      • Count the double bonds along the chain; only accept if exactly two are found.
  If any acid group yields a valid chain, the molecule is accepted.
Note: This is a simplified method and may not capture all edge cases.
"""
from rdkit import Chem

def is_octadecadienoic_acid(smiles: str):
    """
    Determines whether a molecule qualifies as an octadecadienoic acid based on its SMILES string.
    The criteria:
      - Contains at least one carboxylic acid group (matching "C(=O)[OH]" or "O=C(O)").
      - There exists at least one straight, unbranched contiguous chain of 18 carbon atoms (including the acid carbon)
        from the acid carbon to a terminal methyl carbon.
      - All bonds in that chain are either single or double bonds and are not in rings.
      - Exactly 2 of the bonds along the chain are double bonds.
    
    Args:
      smiles (str): SMILES representation of the molecule.
    
    Returns:
      bool: True if the molecule qualifies as an octadecadienoic acid, otherwise False.
      str: Explanation for the decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Step 1. Find carboxylic acid group(s).
    acid_smarts_list = ["C(=O)[OH]", "O=C(O)"]
    acid_matches_all = []
    for smarts in acid_smarts_list:
        pattern = Chem.MolFromSmarts(smarts)
        if pattern is None:
            continue
        matches = mol.GetSubstructMatches(pattern)
        if matches:
            acid_matches_all.extend(matches)
    if not acid_matches_all:
        return False, "No carboxylic acid group found"
    
    # Precompute candidate terminal methyl carbons in the molecule.
    # Terminal methyl: a carbon that has exactly one carbon neighbor.
    candidate_terminals = set()
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:
            carbon_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 6]
            if len(carbon_neighbors) == 1:
                candidate_terminals.add(atom.GetIdx())
    
    if not candidate_terminals:
        return False, "No candidate terminal methyl group found"
    
    target_chain_len = 18  # total number of carbon atoms that must be present in a valid chain
    
    # DFS to look for a chain from start to terminal with exactly target_chain_len carbons.
    def dfs_paths(current, target, max_len, path, paths):
        if len(path) == max_len:
            if current == target:
                paths.append(path.copy())
            return
        # Explore neighbors that are carbons and have not been visited.
        for nbr in mol.GetAtomWithIdx(current).GetNeighbors():
            if nbr.GetAtomicNum() != 6:
                continue
            nbr_idx = nbr.GetIdx()
            if nbr_idx in path:
                continue
            bond = mol.GetBondBetweenAtoms(current, nbr_idx)
            if bond is None:
                continue
            # Only allow single or double bonds and that are not in rings.
            if bond.GetBondType() not in (Chem.BondType.SINGLE, Chem.BondType.DOUBLE):
                continue
            if bond.IsInRing():
                continue
            # Proceed recursively.
            path.append(nbr_idx)
            dfs_paths(nbr_idx, target, max_len, path, paths)
            path.pop()
    
    # For debugging: store details on why paths failed.
    debug_reasons = []
    
    # Try each acid group as the starting point.
    for match in acid_matches_all:
        # Identify the acid carbon within the match (atom with atomic number 6).
        acid_carbon = None
        for idx in match:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() == 6:
                acid_carbon = idx
                break
        if acid_carbon is None:
            debug_reasons.append("Acid group found but no acid carbon identified in one match")
            continue
        
        # For this acid carbon, search among candidate terminal methyl groups.
        for term in candidate_terminals:
            paths = []
            dfs_paths(acid_carbon, term, target_chain_len, [acid_carbon], paths)
            for p in paths:
                # Make sure the path length is exactly 18.
                if len(p) != target_chain_len:
                    continue
                # Check that the chain is linear: each carbon in the chain should not have an extra carbon neighbor
                # that is not part of the chain.
                branched = False
                for atom_idx in p:
                    atom = mol.GetAtomWithIdx(atom_idx)
                    for neighbor in atom.GetNeighbors():
                        if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() not in p:
                            branched = True
                            break
                    if branched:
                        break
                if branched:
                    debug_reasons.append(f"Chain {p} is branched (extra carbon substituent found).")
                    continue
                
                # Count the number of double bonds in consecutive pairs along the chain.
                dbl_bond_count = 0
                valid_chain = True
                for i in range(len(p) - 1):
                    bond = mol.GetBondBetweenAtoms(p[i], p[i+1])
                    if bond is None:
                        valid_chain = False
                        break
                    if bond.GetBondType() == Chem.BondType.DOUBLE:
                        dbl_bond_count += 1
                if not valid_chain:
                    debug_reasons.append(f"Chain {p} has an invalid bond along it.")
                    continue
                if dbl_bond_count != 2:
                    debug_reasons.append(f"Chain {p} has {dbl_bond_count} C=C bond(s), expected exactly 2.")
                    continue
                
                # If we reach here, we have found a valid straight-chain.
                return True, "Molecule is a straight-chain C18 fatty acid with exactly 2 C=C bonds"
    
    if debug_reasons:
        return False, "No valid straight-chain C18 found with exactly 2 C=C bonds. " + "; ".join(debug_reasons)
    else:
        return False, "No valid straight-chain C18 found with exactly 2 C=C bonds."

# Example usage:
if __name__ == "__main__":
    # Test with one provided example: 9-cis,11-trans-octadecadienoic acid
    test_smiles = "CCCCCC\\C=C\\C=C/CCCCCCCC(O)=O"
    result, reason = is_octadecadienoic_acid(test_smiles)
    print(result, reason)