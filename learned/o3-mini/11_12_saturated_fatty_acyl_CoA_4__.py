"""
Classifies: CHEBI:84948 11,12-saturated fatty acyl-CoA(4-)
"""
#!/usr/bin/env python3
"""
Classifies: CHEBI:84947 – 11,12-saturated fatty acyl-CoA(4-)
Definition:
   Any fatty acyl-CoA(4-) in which the 11-12 bond of the fatty acyl group is saturated.
   (The numbering of the fatty acyl group is defined from the carbonyl carbon as C1.)
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_11_12_saturated_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines whether the given SMILES string corresponds to an 11,12-saturated fatty acyl-CoA(4-).
    We first check that the molecule contains a CoA-like fragment and a thioester linkage.
    Then we follow the fatty acyl chain in a strictly linear fashion (rejecting chains with branches)
    starting from the carbonyl carbon (C1). Finally, we check that the bond between C11 and C12 (numbered from C1)
    is a single bond; that is, the 11-12 bond is saturated.
    
    Args:
        smiles (str): SMILES string of the molecule.
    Returns:
        bool: True if the molecule meets the criteria, False otherwise.
        str: An explanation message.
    """
    # Parse SMILES.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for a CoA-like fragment.
    # We use a simple SMARTS pattern to look for a portion of the CoA moiety.
    coa_pattern = Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "CoA moiety not found"
    
    # Look for the thioester group that connects the fatty acyl chain to CoA.
    # SMARTS pattern for a thioester: a carbonyl carbon attached to a sulfur.
    thioester_pattern = Chem.MolFromSmarts("[CX3](=O)[S]")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "Thioester (acyl-CoA linkage) not found"
    
    # Assume the first match is the fatty acyl-CoA connection.
    # In the match, the first atom is the carbonyl carbon (designated as C1).
    thioester_match = thioester_matches[0]
    acyl_carbon_idx = thioester_match[0]
    acyl_carbon = mol.GetAtomWithIdx(acyl_carbon_idx)
    
    # Find the neighbor carbon (ignoring oxygen and sulfur atoms) that starts the fatty acyl chain.
    acyl_chain_start = None
    for nbr in acyl_carbon.GetNeighbors():
        if nbr.GetAtomicNum() == 6:
            acyl_chain_start = nbr
            break
    if acyl_chain_start is None:
        return False, "Fatty acyl chain not found (no carbon attached to acyl carbon)"
    
    # Follow a linear (unbranched) carbon chain starting from acyl_chain_start.
    # We construct a list 'chain' starting with the acyl carbon (C1) followed by sequential carbons.
    chain = [acyl_carbon, acyl_chain_start]
    current = acyl_chain_start
    previous = acyl_carbon
    while True:
        # Get all carbon neighbors of the current atom excluding the atom we just came from.
        next_carbons = [nbr for nbr in current.GetNeighbors() if nbr.GetAtomicNum() == 6 and nbr.GetIdx() != previous.GetIdx()]
        if len(next_carbons) == 1:
            chain.append(next_carbons[0])
            previous, current = current, next_carbons[0]
        else:
            # No extension or branching encountered: stop following the chain.
            break

    # To ensure we are dealing with a straight-chain fatty acyl group,
    # check that each intermediate carbon (except the first and last) has no extra carbon branch.
    for i in range(1, len(chain) - 1):
        atom = chain[i]
        # The two chain neighbors are already part of the main chain.
        chain_neighbors = {chain[i - 1].GetIdx(), chain[i + 1].GetIdx()}
        extra_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 6 and nbr.GetIdx() not in chain_neighbors]
        if extra_neighbors:
            return False, "Fatty acyl chain appears branched"
    
    # Ensure the chain is long enough for an 11-12 bond check
    if len(chain) < 12:
        return False, f"Fatty acyl chain too short (length {len(chain)}); at least 12 carbons required for 11-12 bond check."
    
    # By definition, the carbonyl carbon is C1 so the bond between the 11th and 12th carbon
    # is between chain[10] and chain[11] (using 0-indexing).
    bond = mol.GetBondBetweenAtoms(chain[10].GetIdx(), chain[11].GetIdx())
    if bond is None:
        return False, "Bond between C11 and C12 not found"
    if bond.GetBondType() != Chem.BondType.SINGLE:
        return False, "11-12 bond is not a single (saturated) bond"
    
    return True, "Molecule is an 11,12-saturated fatty acyl-CoA(4-) with a proper fatty acyl chain and CoA moiety"


# For testing (the following block may be removed in production)
if __name__ == "__main__":
    # Example test: (3R,17Z,20Z,23Z,26Z)-3-hydroxydotriacontatetraenoyl-CoA(4-)
    test_smiles = r"CCCC\C=C/C\C=C/C\C=C/C\C=C/CCCCCCCCCCCCC[C@@H](O)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)n1cnc2c(N)ncnc12"
    result, reason = is_11_12_saturated_fatty_acyl_CoA_4__(test_smiles)
    print(result, reason)