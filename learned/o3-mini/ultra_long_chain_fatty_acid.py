"""
Classifies: CHEBI:143004 ultra-long-chain fatty acid
"""
"""
Classifies: Ultra Long‐Chain Fatty Acid
Definition: A fatty acid whose free carboxylic acid group (C(=O)[OH]) is attached to a continuous acyclic carbon chain 
            (possibly with minor branching) that has a length (counting the acid carbon) greater than 27.
            In addition the molecule must be “pure” (i.e. nearly all non‐H atoms belong to the fatty acid fragment).
            
The approach:
  1. Look for a free carboxylic acid group via SMARTS "[CX3](=O)[OX2H1]".
  2. Ensure that the acid carbon has exactly one carbon neighbor (the fatty acyl chain start).
  3. Using a DFS (ignoring atoms in rings), follow only carbon neighbors to obtain the longest connected acyclic carbon set.
  4. Form the “fatty acid fragment” as that set plus the acid carbon and its two oxygen atoms.
  5. Compare the number of heavy atoms in this fragment to the molecule’s total heavy atoms, and only accept
     if most (≥80%) belong to the fragment – this reduces false positives from complex molecules.
  6. Finally, if the computed chain length (fatty acyl chain plus the acid carbon) is > 27, then classify it as an ultra‐long‐chain fatty acid.
"""

from rdkit import Chem

def is_ultra_long_chain_fatty_acid(smiles: str):
    """
    Determines whether a molecule qualifies as an ultra long-chain fatty acid.
    A qualifying molecule must have a free carboxylic acid group (C(=O)[OH])
    where the acid carbon is attached to a continuous, acyclic chain of carbons such that 
    (chain length counting the acid carbon) > 27. Also the molecule should be largely
    composed of the fatty acyl fragment so that we do not falsely classify complex molecules.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule qualifies, False otherwise.
        str: Explanation of the classification.
    """
    # Parse input SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for free carboxylic acid group: C(=O)[OH]
    acid_smarts = "[CX3](=O)[OX2H1]"
    acid_pattern = Chem.MolFromSmarts(acid_smarts)
    acid_matches = mol.GetSubstructMatches(acid_pattern)
    if not acid_matches:
        return False, "No free carboxylic acid group found"
    
    # We assume that the molecule should behave as a fatty acid if one (or more) acid groups are present.
    # We will try each acid group and classify if one qualifies.
    fatty_acid_found = False
    explanations = []
    
    # Helper DFS: Given a carbon atom, traverse neighbors that are carbon, not in a ring.
    def dfs_carbon(atom, visited):
        visited.add(atom.GetIdx())
        length = 1  # count current carbon
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6 and nbr.GetIdx() not in visited and not nbr.IsInRing():
                # continue DFS on this carbon neighbor
                branch_len = 1 + dfs_carbon(nbr, visited.copy())
                if branch_len > length:
                    length = branch_len
        return length

    # We also want to collect the set of atom indices that are part of the acyclic carbon chain.
    def dfs_collect(atom, visited):
        visited.add(atom.GetIdx())
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6 and nbr.GetIdx() not in visited and not nbr.IsInRing():
                dfs_collect(nbr, visited)
        return visited

    # Count heavy atoms in the molecule (atoms with atomic number > 1)
    total_heavy = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1)
    
    # Process each acid match
    for match in acid_matches:
        # According to our SMARTS "[CX3](=O)[OX2H1]", index 0 is the acid (C) and index 1 is the hydroxyl oxygen.
        acid_c = mol.GetAtomWithIdx(match[0])
        
        # The acid carbon should have exactly one carbon neighbor (the fatty acyl chain start)
        chain_neighs = [nbr for nbr in acid_c.GetNeighbors() if nbr.GetAtomicNum() == 6]
        if len(chain_neighs) != 1:
            explanations.append(f"Acid carbon (idx {acid_c.GetIdx()}) does not have exactly one carbon neighbor, found {len(chain_neighs)}")
            continue
        
        chain_start = chain_neighs[0]
        # Calculate longest acyclic carbon chain starting from chain_start.
        chain_length = dfs_carbon(chain_start, set())
        # Total chain length including the acid carbon itself.
        total_chain_atoms = chain_length + 1
        
        # Collect all acyclic carbons in the fatty acyl fragment (the DFS may be branched)
        chain_atom_idxs = dfs_collect(chain_start, set())
        # Build the fatty acid fragment: include the acid carbon and the two oxygens directly bound to it.
        fragment_idxs = set(chain_atom_idxs)
        fragment_idxs.add(acid_c.GetIdx())
        # Add oxygens attached to acid carbon that form the acid group.
        for nbr in acid_c.GetNeighbors():
            if nbr.GetAtomicNum() == 8:
                fragment_idxs.add(nbr.GetIdx())
        
        # Count heavy atoms in the fatty acid fragment.
        frag_heavy = sum(1 for idx in fragment_idxs if mol.GetAtomWithIdx(idx).GetAtomicNum() > 1)
        
        # We require that the fatty acid fragment makes up most (>=80%) of the heavy atoms.
        ratio = frag_heavy / total_heavy
        if ratio < 0.8:
            explanations.append(f"Fatty acyl fragment only accounts for {ratio:.0%} of heavy atoms; molecule appears to have substantial extra substructures")
            continue
        
        # Now, if overall chain length exceeds 27 carbons, we classify as ultra-long-chain.
        if total_chain_atoms > 27:
            return True, f"Chain length is {total_chain_atoms} carbons, which qualifies as an ultra-long-chain fatty acid; fatty acid fragment comprises {ratio:.0%} of heavy atoms."
        else:
            explanations.append(f"Chain length is {total_chain_atoms} carbons, which does not exceed the C27 threshold")
    
    # If none of the acid groups qualified, return False.
    if explanations:
        return False, " ; ".join(explanations)
    return False, "No qualifying free carboxylic acid with a sufficiently long and pure chain was found"

# Example usage:
if __name__ == "__main__":
    # Test with dotriacontanoic acid (32-carbon fatty acid)
    test_smiles = "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC(O)=O"
    result, explanation = is_ultra_long_chain_fatty_acid(test_smiles)
    print(result, explanation)