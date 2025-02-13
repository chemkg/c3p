"""
Classifies: CHEBI:33184 long-chain fatty acyl-CoA
"""
#!/usr/bin/env python3
"""
Classifies: long-chain fatty acyl-CoA 
(a fatty acyl-CoA derived from the condensation of a long-chain (C13 to C22) fatty acid with coenzyme A)
"""

from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

def is_long_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a long-chain fatty acyl-CoA based on its SMILES string.
    
    It does so by checking for a CoA moiety (using an adenine-based fragment as marker)
    and a thioester group (C(=O)S) linking a fatty acyl chain. The fatty acyl chain, including
    the carbonyl carbon, must contain between 13 and 22 carbon atoms.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a long-chain fatty acyl-CoA, False otherwise
        str: Reason for the classification
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a SMARTS pattern for a CoA moiety.
    # We search for an adenine fragment present in coenzyme A.
    # Use a more complete adenine SMARTS pattern.
    coa_substructure_smarts = "c1nc2c(n1)nc(n2)"  # adenine fragment
    coa_frag = Chem.MolFromSmarts(coa_substructure_smarts)
    if coa_frag is None:
        return False, "Invalid CoA SMARTS pattern"
    if not mol.HasSubstructMatch(coa_frag):
        return False, "No CoA (coenzyme A) moiety detected"
    
    # Define a SMARTS pattern for a thioester group (C(=O)S).
    thioester_smarts = "[CX3](=O)[SX2]"
    thioester = Chem.MolFromSmarts(thioester_smarts)
    if thioester is None:
        return False, "Invalid thioester SMARTS pattern"
    thioester_matches = mol.GetSubstructMatches(thioester)
    if not thioester_matches:
        return False, "No thioester group (linking fatty acid to CoA) detected"

    # For each thioester match, identify and traverse the potential fatty acyl chain.
    valid_match_found = False
    reasons = []
    for match in thioester_matches:
        # The matching returns indices for (carbonyl carbon, sulfur)
        carbonyl_idx, sulfur_idx = match
        carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
        # Out of the neighbors of the carbonyl carbon, one is the double-bonded oxygen and one is sulfur.
        # The remaining neighbor (usually a carbon) is the start of the acyl chain.
        acyl_neighbors = []
        for neighbor in carbonyl_atom.GetNeighbors():
            # Exclude oxygen (likely the double-bonded oxygen) and the sulfur
            if neighbor.GetIdx() == sulfur_idx:
                continue
            if neighbor.GetAtomicNum() == 6:  # only consider carbon atoms
                acyl_neighbors.append(neighbor)
        if len(acyl_neighbors) != 1:
            reasons.append("Ambiguous or missing acyl chain for one thioester match")
            continue
        acyl_start = acyl_neighbors[0]

        # Perform a depth-first search (restricted to C-C bonds) to collect all carbon atoms
        # that are part of the fatty acyl chain.
        # Include the carbonyl carbon in the chain count.
        acyl_atoms = set()
        acyl_atoms.add(carbonyl_idx)
        stack = [acyl_start.GetIdx()]
        while stack:
            curr_idx = stack.pop()
            if curr_idx in acyl_atoms:
                continue
            curr_atom = mol.GetAtomWithIdx(curr_idx)
            if curr_atom.GetAtomicNum() != 6:
                continue
            acyl_atoms.add(curr_idx)
            for nb in curr_atom.GetNeighbors():
                # Only traverse bonds between carbons.
                if nb.GetAtomicNum() == 6 and nb.GetIdx() not in acyl_atoms:
                    stack.append(nb.GetIdx())
        
        fatty_acyl_chain_length = len(acyl_atoms)
        # Check if the fatty acyl chain length (including the carbonyl carbon) is between 13 and 22.
        if 13 <= fatty_acyl_chain_length <= 22:
            valid_match_found = True
            return True, (f"Found thioester linking a fatty acyl chain with {fatty_acyl_chain_length} carbons "
                          f"(within expected range 13-22) attached to a CoA moiety")
        else:
            reasons.append(f"Fatty acyl chain with {fatty_acyl_chain_length} carbons; expected between 13 and 22")

    # If we did not find any valid thioester+acyl chain match then return failure with reasons.
    if not valid_match_found:
        reason_str = "; ".join(reasons) if reasons else "No appropriate thioester fatty acyl chain found"
        return False, reason_str

    return False, "No valid long-chain fatty acyl-CoA structure found"

# This module can be further enhanced by unit testing with the provided SMILES examples.