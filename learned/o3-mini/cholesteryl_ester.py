"""
Classifies: CHEBI:17002 cholesteryl ester
"""
"""
Classifies: cholesteryl ester
Definition: A sterol ester obtained by formal condensation of the carboxy group of any carboxylic acid 
with the 3-hydroxy group of cholesterol.
Heuristic: The molecule must contain one ester group (O–C(=O)) whose alkoxy oxygen is attached to 
a fused ring system (at least 4 rings, typical for a steroid nucleus) and whose carbonyl carbon is attached 
to a non-cyclic (fatty acid) chain of sufficient length.
"""

from rdkit import Chem

def is_cholesteryl_ester(smiles: str):
    """
    Determines if a molecule is a cholesteryl ester based on its SMILES string.
    
    The algorithm checks:
      1. The SMILES is valid.
      2. The molecule contains exactly one ester group of the form O–C(=O).
      3. The oxygen (alkoxy group) is attached to a fused ring system (at least 4 rings)
         that is typical of the cholesterol (steroid) moiety.
      4. The carbonyl carbon (C=O) is attached to a non-cyclic, aliphatic chain
         (representing the fatty acid part) with a minimum chain length.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule is classified as a cholesteryl ester, False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a SMARTS query for an ester fragment.
    # This query looks for an alkoxy oxygen (with no hydrogen) attached to a carbonyl carbon.
    ester_query = Chem.MolFromSmarts("[O;!H0][C](=O)")
    ester_matches = mol.GetSubstructMatches(ester_query)
    if len(ester_matches) != 1:
        return False, "Molecule must contain exactly one ester group"

    # Get the matched atoms: by convention, match[0] is the ester oxygen, match[1] the carbonyl carbon.
    ester_match = ester_matches[0]
    oxy_idx = ester_match[0]
    carbonyl_idx = ester_match[1]
    atom_oxy = mol.GetAtomWithIdx(oxy_idx)
    atom_caro = mol.GetAtomWithIdx(carbonyl_idx)

    # Find the neighbor of the ester oxygen that is not the carbonyl carbon.
    # This should correspond to the cholesterol part (from its 3-hydroxy group).
    oxy_neighbors = atom_oxy.GetNeighbors()
    chol_neighbor = None
    for nb in oxy_neighbors:
        if nb.GetIdx() != carbonyl_idx:
            chol_neighbor = nb
            break
    if chol_neighbor is None:
        return False, "Ester oxygen is not attached to any substituent (cholesterol part missing)"
        
    # Verify that the molecule contains a fused ring system.
    # Cholesterol has a tetracyclic (four-ring) steroid nucleus.
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    if len(rings) < 4:
        return False, "The molecule does not contain enough rings (expected count for a steroid nucleus is 4)"
    
    # Check that the group attached via the ester oxygen is indeed part of a ring system.
    if not chol_neighbor.IsInRing():
        return False, "The substituent attached to the ester oxygen is not part of a fused ring (steroid) system"

    # Next, verify the fatty acyl part.
    # In an ester, the carbonyl carbon is attached to one oxygen (already used) 
    # and another neighbor that should be the fatty acid chain.
    fatty_neighbor = None
    for nb in atom_caro.GetNeighbors():
        if nb.GetIdx() != atom_oxy.GetIdx():
            fatty_neighbor = nb
            break
    if fatty_neighbor is None:
        return False, "Carbonyl carbon does not have a fatty acid chain attached"

    # We expect the fatty acid chain to be acyclic. If this neighbor is in a ring, it is unexpected.
    if fatty_neighbor.IsInRing():
        return False, "The fatty acid chain appears to be cyclic rather than a linear aliphatic chain"

    # Define a helper function to perform a depth-first search starting from the fatty acid attachment.
    # This function walks only through non-ring, carbon atoms and returns the maximum chain length (in bonds).
    def get_max_chain_length(atom, visited):
        max_len = 0
        for nbr in atom.GetNeighbors():
            if nbr.GetIdx() in visited:
                continue
            if nbr.IsInRing():  # avoid stepping into cyclic parts
                continue
            if nbr.GetAtomicNum() != 6:  # consider only carbon atoms
                continue
            new_visited = visited | {nbr.GetIdx()}
            length = 1 + get_max_chain_length(nbr, new_visited)
            if length > max_len:
                max_len = length
        return max_len

    # Compute the chain length from the fatty acid attachment.
    chain_length = get_max_chain_length(fatty_neighbor, {fatty_neighbor.GetIdx()})
    # We require at least a chain of 5 bonds (i.e. at least 6 carbon atoms in a row) to be considered fatty.
    if chain_length < 5:
        return False, f"Fatty acid chain too short (chain length = {chain_length + 1} carbons expected ≥ 6)"

    return True, "Molecule is a cholesteryl ester with a steroid nucleus and a fatty acid chain attached via an ester linkage"

# Example usage:
if __name__ == "__main__":
    # Test one of the provided examples:
    test_smiles = "O([C@@H]1CC=2[C@@](C3C(C4[C@@](C(CC4)[C@@H](CCCC(C)C)C)(CC3)C)CC2)(CC1)C)C(CCCCCCCCCCC/C=C\\CCCCCCCC)=O"
    result, reason = is_cholesteryl_ester(test_smiles)
    print(result, reason)