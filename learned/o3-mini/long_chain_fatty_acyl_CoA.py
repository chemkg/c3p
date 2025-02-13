"""
Classifies: CHEBI:33184 long-chain fatty acyl-CoA
"""
#!/usr/bin/env python3
"""
Classifies: long-chain fatty acyl-CoA (a fatty acyl-CoA derived from condensation of a long-chain (C13 to C22) fatty acid with coenzyme A)
"""

from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

def is_long_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a long-chain fatty acyl-CoA based on its SMILES string.
    It does so by checking for a CoA moiety (using an adenine fragment as marker)
    and a thioester (C(=O)S) linking a fatty acyl chain. The fatty acyl chain (including the carbonyl carbon)
    must contain between 13 and 22 carbon atoms.

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

    # Check for a CoA moiety.
    # We look for an adenine-like substructure present in Coenzyme A.
    coa_substructure_smarts = "n1cnc2"  # adenine fragment
    coa_frag = Chem.MolFromSmarts(coa_substructure_smarts)
    if not mol.HasSubstructMatch(coa_frag):
        return False, "No CoA (coenzyme A) moiety detected"

    # Search for a thioester group (carbonyl attached to sulfur).
    thioester_smarts = "[CX3](=O)[SX2]"
    thioester = Chem.MolFromSmarts(thioester_smarts)
    thioester_matches = mol.GetSubstructMatches(thioester)
    if not thioester_matches:
        return False, "No thioester group (linking fatty acid to CoA) detected"

    # For each thioester match look for an acyl chain with the appropriate number of carbons.
    # The thioester match returns a tuple of atom indices: (carbonyl_carbon, sulfur)
    for match in thioester_matches:
        carbonyl_idx, sulfur_idx = match
        carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
        # Identify the acyl chain neighbor.
        # The carbonyl carbon normally has three neighbors:
        #   1) one double-bonded oxygen (should be ignored),
        #   2) the sulfur atom (link to CoA, ignore),
        #   3) the fatty acyl carbon.
        acyl_neighbors = []
        for neighbor in carbonyl_atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 6:  # only consider carbon atoms
                acyl_neighbors.append(neighbor)
        if len(acyl_neighbors) != 1:
            # Either ambiguous or missing acyl chain; try the next match.
            continue
        acyl_start = acyl_neighbors[0]

        # Now perform a DFS (restricted to C-C bonds) to collect all carbon atoms
        # that are part of the fatty acyl chain.
        # We include the carbonyl carbon as part of the fatty acid.
        acyl_atoms = set()
        acyl_atoms.add(carbonyl_idx)
        stack = [acyl_start.GetIdx()]
        while stack:
            curr_idx = stack.pop()
            if curr_idx in acyl_atoms:
                continue
            curr_atom = mol.GetAtomWithIdx(curr_idx)
            # Only traverse along carbon atoms
            if curr_atom.GetAtomicNum() != 6:
                continue
            acyl_atoms.add(curr_idx)
            for nb in curr_atom.GetNeighbors():
                if nb.GetAtomicNum() == 6 and nb.GetIdx() not in acyl_atoms:
                    # Do not cross into the CoA moiety: the thioester sulfur's branch is not followed
                    stack.append(nb.GetIdx())

        fatty_acyl_chain_length = len(acyl_atoms)
        # The fatty acyl chain includes the carbonyl carbon.
        if 13 <= fatty_acyl_chain_length <= 22:
            return True, (f"Found thioester linking a fatty acyl chain with {fatty_acyl_chain_length} carbons "
                          f"(within the expected range of 13-22) attached to a CoA moiety")
        else:
            return False, (f"Fatty acyl chain has {fatty_acyl_chain_length} carbons; expected between 13 and 22")

    # If no thioester match produced a valid fatty acyl chain, return a negative result.
    return False, "No appropriate thioester fatty acyl chain found"
    
# Note: This module can be further enhanced by testing with the provided SMILES examples.