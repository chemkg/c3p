"""
Classifies: CHEBI:17002 cholesteryl ester
"""
"""
Classifies: cholesteryl ester
Definition: A sterol ester obtained by formal condensation of the carboxy group of any carboxylic acid 
with the 3-hydroxy group of cholesterol.
Heuristic: The molecule must contain one ester linkage (–O–C(=O)–)
where the alkoxy oxygen is connected to a fused ring system (steroid nucleus, typically ≥4 rings)
and the carbonyl carbon is connected to a non-cyclic, sufficiently long fatty acid chain.
"""

from rdkit import Chem

def is_cholesteryl_ester(smiles: str):
    """
    Determines if a molecule is a cholesteryl ester based on its SMILES string.
    
    The algorithm checks:
      1. That the SMILES string is valid.
      2. That the molecule contains at least one ester group of the form O–C(=O) and then 
         inspects each candidate to see if it fulfills:
         a) The oxygen (alkoxy) is attached to a substituent that is in a ring system and
            the overall molecule has at least 4 rings (typical for the steroid nucleus).
         b) The carbonyl carbon is attached to a fatty acid chain that is acyclic (i.e. not in any ring)
            and has a minimum chain length (≥6 carbons in a row).
      3. If exactly one candidate ester linkage passes these tests, the molecule is classified as a cholesteryl ester.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a cholesteryl ester, False otherwise.
        str: Explanation for the decision.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # First, require that the overall molecule has at least 4 rings (steroid nucleus)
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    if len(rings) < 4:
        return False, "Molecule does not contain enough rings for a steroid nucleus (expected ≥4 rings)"
    
    # Define a SMARTS pattern for an ester group:
    # The pattern looks for an oxygen (without hydrogen) bonded to a carbonyl carbon.
    ester_query = Chem.MolFromSmarts("[O;!H0][C](=O)")
    ester_matches = mol.GetSubstructMatches(ester_query)
    if not ester_matches:
        return False, "No ester group found"
    
    valid_candidates = []
    
    # Helper: compute maximum acyclic carbon chain length starting from a given atom.
    def get_max_chain_length(atom, visited):
        max_len = 0
        for nbr in atom.GetNeighbors():
            if nbr.GetIdx() in visited:
                continue
            # only consider non-ring carbon atoms
            if nbr.IsInRing() or nbr.GetAtomicNum() != 6:
                continue
            new_visited = visited | {nbr.GetIdx()}
            length = 1 + get_max_chain_length(nbr, new_visited)
            if length > max_len:
                max_len = length
        return max_len

    # Iterate over each ester candidate
    for match in ester_matches:
        # By our SMARTS, match[0] is the ester oxygen and match[1] is the carbonyl carbon.
        oxy_idx, carbonyl_idx = match[0], match[1]
        atom_oxy = mol.GetAtomWithIdx(oxy_idx)
        atom_caro = mol.GetAtomWithIdx(carbonyl_idx)
        
        # For the ester oxygen, find the neighbor which is not the carbonyl carbon.
        oxy_neighbors = [nbr for nbr in atom_oxy.GetNeighbors() if nbr.GetIdx() != carbonyl_idx]
        if not oxy_neighbors:
            continue  # No substituent attached to ester oxygen
        chol_candidate = oxy_neighbors[0]
        if not chol_candidate.IsInRing():
            continue  # The substituent should belong to the steroid (ring) system
        
        # For additional assurance, check that at least one ring in the molecule contains the cholesterol part.
        in_ring = False
        for ring in rings:
            if chol_candidate.GetIdx() in ring:
                in_ring = True
                break
        if not in_ring:
            continue
        
        # Now, for the carbonyl carbon, find the fatty acid chain substituent:
        fatty_neighbors = [nbr for nbr in atom_caro.GetNeighbors() if nbr.GetIdx() != oxy_idx]
        if not fatty_neighbors:
            continue
        fatty_candidate = fatty_neighbors[0]
        # The fatty acid branch should be acyclic.
        if fatty_candidate.IsInRing():
            continue
        
        # Determine the chain length of the fatty acid part by walking along connected carbon atoms (non-ring).
        chain_length = get_max_chain_length(fatty_candidate, {fatty_candidate.GetIdx()})
        # Require at least 5 bonds (which means at least 6 carbon atoms in a carbon chain).
        if chain_length < 5:
            continue
        
        # This candidate passes all tests
        valid_candidates.append(match)
        
    if len(valid_candidates) == 0:
        return False, "No valid cholesteryl ester linkage found: ester group missing or not in proper cholesterol/fatty acid context"
    elif len(valid_candidates) > 1:
        return False, f"Multiple ({len(valid_candidates)}) ester linkages found; expected exactly one cholesteryl ester linkage"
    else:
        return True, "Molecule is classified as a cholesteryl ester with a steroid nucleus and a valid fatty acid chain"

# Example usage:
if __name__ == "__main__":
    # Test one example from the provided list:
    test_smiles = "O([C@@H]1CC=2[C@@](C3C(C4[C@@](C(CC4)[C@@H](CCCC(C)C)C)(CC3)C)CC2)(CC1)C)C(CCCCCCCCCCC/C=C\\CCCCCCCC)=O"
    result, reason = is_cholesteryl_ester(test_smiles)
    print(result, reason)