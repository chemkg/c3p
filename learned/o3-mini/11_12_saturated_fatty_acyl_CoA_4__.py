"""
Classifies: CHEBI:84948 11,12-saturated fatty acyl-CoA(4-)
"""
#!/usr/bin/env python3
"""
Classifies: CHEBI:84947 partner – 11,12-saturated fatty acyl-CoA(4-)
Definition:
   Any fatty acyl-CoA(4-) in which the 11-12 bond of the fatty acyl group is saturated.
   (The numbering of the fatty acyl group is defined from the carbonyl carbon as C1.)
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_11_12_saturated_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines whether the given SMILES string corresponds to an 11,12-saturated fatty acyl-CoA(4-).
    This is defined as a fatty acyl-CoA (i.e. a molecule with a thioester linkage attached to a CoA-like fragment)
    for which the fatty acyl chain (starting at the carbonyl carbon) has its 11-12 bond (i.e. between the 11th and 12th carbon)
    saturated (single bond).

    Args:
        smiles (str): SMILES string of the molecule.
    Returns:
        bool: True if the molecule is an 11,12-saturated fatty acyl-CoA(4-), else False.
        str: A message explaining the basis of the classification.
    """
    # Parse SMILES:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for a CoA fragment.
    # This is a crude substructure search: many acyl-CoA molecules contain the motif "SCCNC(=O)CCNC(=O)"
    coa_pattern = Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "CoA moiety not found"

    # Look for the thioester group that connects the fatty acyl chain to CoA.
    # We use the SMARTS pattern for a thioester: a carbonyl carbon attached to a sulfur.
    thioester_pattern = Chem.MolFromSmarts("[CX3](=O)[S]")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "Thioester (acyl-CoA linkage) not found"

    # Assume the first thioester match is our fatty acyl-CoA connection.
    # In the pattern "[CX3](=O)[S]" the first atom is the carbonyl carbon.
    thioester_match = thioester_matches[0]
    acyl_carbon = thioester_match[0]  # Carbonyl carbon (C1 of the acyl chain)
    
    # Identify the neighbor of the acyl carbon that belongs to the fatty acyl chain.
    # The carbonyl carbon in a fatty acid is bonded to the carbonyl oxygen and to the sulfur.
    # Its additional neighbor (if any) should be the alpha carbon of the fatty acyl chain.
    atom_acyl = mol.GetAtomWithIdx(acyl_carbon)
    acyl_chain_start = None
    for neighbor in atom_acyl.GetNeighbors():
        # Skip oxygen and sulfur atoms
        if neighbor.GetAtomicNum() == 6:
            acyl_chain_start = neighbor
            break
    if acyl_chain_start is None:
        return False, "Fatty acyl chain not found (no carbon attached to acyl carbon)"

    # Now we have to follow the main (longest) carbon chain that constitutes the fatty acyl chain.
    # We define a recursive DFS that only walks through carbon atoms, ignoring rings and branches that
    # might represent side chains.
    def dfs(atom, visited):
        paths = []
        visited.add(atom.GetIdx())
        extended = False
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() != 6:
                continue  # only follow carbon atoms
            if nbr.GetIdx() in visited:
                continue
            # Call recursively
            subpaths = dfs(nbr, visited.copy())
            for sp in subpaths:
                paths.append([atom] + sp)
            extended = True
        if not extended:
            # no further extension: return a path with only this atom
            return [[atom]]
        return paths

    # We start our DFS from the acyl_chain_start and prepend the acyl carbon to have the complete chain.
    all_paths = dfs(acyl_chain_start, set())
    if not all_paths:
        return False, "No carbon chain found in fatty acyl group"
    # Prepend the acyl carbon (C1) to each path:
    all_paths = [[mol.GetAtomWithIdx(acyl_carbon)] + path for path in all_paths]
    # Choose the longest path (by number of atoms)
    longest_path = max(all_paths, key=lambda p: len(p))
    chain_length = len(longest_path)
    if chain_length < 12:
        return False, f"Fatty acyl chain too short (length {chain_length}). At least 12 carbons required for 11-12 bond check."

    # Check the 11-12 bond: numbering from the carbonyl carbon (which is longest_path[0])
    # means that the bond between atoms at positions 10 and 11 (0-indexed) is the 11-12 bond.
    idx1 = longest_path[10].GetIdx()
    idx2 = longest_path[11].GetIdx()
    bond = mol.GetBondBetweenAtoms(idx1, idx2)
    if bond is None:
        return False, "Bond between C11 and C12 not found in the chain"
    # The bond should be a single bond to be considered saturated.
    if bond.GetBondType() != Chem.BondType.SINGLE:
        return False, "11-12 bond is not a single (saturated) bond"
    
    return True, "Molecule is an 11,12-saturated fatty acyl-CoA(4-) with a proper fatty acyl chain and CoA moiety"

# For testing (you can remove or comment out this section in production)
if __name__ == "__main__":
    # Example: (3R,17Z,20Z,23Z,26Z)-3-hydroxydotriacontatetraenoyl-CoA(4-)
    test_smiles = r"CCCC\C=C/C\C=C/C\C=C/C\C=C/CCCCCCCCCCCCC[C@@H](O)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)n1cnc2c(N)ncnc12"
    result, reason = is_11_12_saturated_fatty_acyl_CoA_4__(test_smiles)
    print(result, reason)