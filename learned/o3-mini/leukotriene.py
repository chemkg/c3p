"""
Classifies: CHEBI:25029 leukotriene
"""
#!/usr/bin/env python
"""
Classifies: leukotriene
Definition:
  Any icosanoid derived from arachidonic acid (C20 polyunsaturated fatty acid) and its derivatives,
  with a backbone that has four carbon–carbon double bonds (i.e. unsaturations) in which three appear as
  a contiguous, alternating (conjugated) sequence.
  
This program attempts to first identify the longest continuous carbon chain in the molecule.
It then checks that the chain is at least 20 carbons long, has exactly 4 C–C double bonds (along that chain),
and contains at least one contiguous segment of three conjugated double bonds in the pattern: C=C-single-C=C-single-C=C.
"""

from rdkit import Chem

def longest_carbon_chain(mol):
    """
    Given an RDKit molecule, returns the longest simple path (list of atom indices)
    that passes only through carbon atoms.
    """
    # Get the indices of all carbon atoms.
    carbon_indices = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    # Build a dictionary mapping carbon atom index to its connected carbon neighbors.
    carbon_neighbors = {}
    for idx in carbon_indices:
        atom = mol.GetAtomWithIdx(idx)
        neighbors = []
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6:
                neighbors.append(nbr.GetIdx())
        carbon_neighbors[idx] = neighbors

    # Use DFS to search for the longest simple path
    longest_path = []
    
    def dfs(current, path, visited):
        nonlocal longest_path
        # Update longest path if needed.
        if len(path) > len(longest_path):
            longest_path = path[:]
        for nbr in carbon_neighbors.get(current, []):
            if nbr not in visited:
                visited.add(nbr)
                dfs(nbr, path + [nbr], visited)
                visited.remove(nbr)
                
    for start in carbon_indices:
        dfs(start, [start], set([start]))
    
    return longest_path

def is_leukotriene(smiles: str):
    """
    Determines if a molecule is a leukotriene based on its SMILES string.
    The criteria used here are:
      1. The molecule must be parseable.
      2. Its longest continuous carbon chain (backbone) must have at least 20 carbons.
      3. Along that chain there must be exactly 4 carbon–carbon double bonds.
      4. Along the chain, there must be at least one contiguous segment in which three
         alternating (conjugated) double bonds appear, i.e. a pattern of bonds: D–S–D–S–D,
         where D is a double and S is a single bond.
         
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule meets the leukotriene criteria, False otherwise.
        str: Explanation of the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Step 1. Identify the longest continuous carbon chain.
    chain = longest_carbon_chain(mol)
    if len(chain) < 20:
        return False, f"Longest carbon chain has {len(chain)} carbons; expected at least 20 for an icosanoid backbone."

    # Step 2. Count the C–C double bonds along the chain.
    double_bond_count = 0
    chain_bond_types = []  # record bond types along the chain in order (for later pattern search)
    for i in range(len(chain)-1):
        bond = mol.GetBondBetweenAtoms(chain[i], chain[i+1])
        if bond is None:
            # This should not happen because chain was built from connected carbons
            chain_bond_types.append("NA")
            continue
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
            double_bond_count += 1
            chain_bond_types.append("D")
        else:
            chain_bond_types.append("S")
    
    if double_bond_count != 4:
        return False, f"Found {double_bond_count} carbon–carbon double bonds along the longest chain; expected exactly 4."

    # Step 3. Look for a contiguous segment of three conjugated double bonds.
    # That means we need a sequence of five bonds in the chain that follow the pattern: D, S, D, S, D.
    found_conjugated = False
    # There are len(chain_bond_types) bonds; sliding window of size 5.
    for i in range(len(chain_bond_types) - 4):
        segment = chain_bond_types[i:i+5]
        if segment == ["D","S","D","S","D"]:
            found_conjugated = True
            break
    if not found_conjugated:
        return False, "No contiguous conjugated tri-double-bond (pattern D-S-D-S-D) system found along the longest carbon chain."
        
    return True, ("Molecule has a longest carbon backbone of {} carbons, with exactly 4 C=C bonds along it, "
                  "including a contiguous conjugated system of three double bonds typical of leukotrienes."
                 ).format(len(chain))

# Example test cases (uncomment to test)
# test_smiles = "CCCC\\C=C/C[C@@H](O)\\C=C\\C=C\\C=C\\[C@@H](O)CCCC(O)=O"  # 6-trans-leukotriene B4
# print(is_leukotriene(test_smiles))