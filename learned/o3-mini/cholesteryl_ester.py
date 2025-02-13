"""
Classifies: CHEBI:17002 cholesteryl ester
"""
"""
Classifies: cholesteryl ester
Definition: A sterol ester obtained by formal condensation of the carboxy group of any carboxylic acid 
with the 3-hydroxy group of cholesterol.
Heuristic: The molecule must contain one ester linkage (–O–C(=O)–)
where the alkoxy oxygen is connected to a fused ring system (steroid nucleus, typically ≥4 rings)
and the carbonyl carbon is connected to an acyclic fatty acid chain of sufficient length.
"""

from rdkit import Chem

def is_cholesteryl_ester(smiles: str):
    """
    Determines if a molecule is a cholesteryl ester based on its SMILES string.
    
    The algorithm checks:
      1. Validity of the SMILES string.
      2. That the overall molecule contains a steroid nucleus (≥4 rings).
      3. That the molecule has at least one ester group matching the pattern O[C](=O).
      4. For each ester candidate:
         a. The oxygen (alkoxy) must be attached to a substituent that is in a ring 
            (consistent with the cholesterol nucleus).
         b. The carbonyl carbon must be attached to an acyclic fatty acid chain
            that is sufficiently long (≥5 successive bonds connecting non‐ring carbons).
      5. Exactly one candidate ester linkage must pass these tests.
    
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
    
    # Require that the overall molecule has at least 4 rings (sterol nucleus)
    rings = mol.GetRingInfo().AtomRings()
    if len(rings) < 4:
        return False, "Molecule does not contain enough rings for a steroid nucleus (expected ≥4 rings)"
    
    # Define SMARTS for an ester linkage: O[C](=O)
    ester_query = Chem.MolFromSmarts("O[C](=O)")
    ester_matches = mol.GetSubstructMatches(ester_query)
    if not ester_matches:
        return False, "No ester group found"
    
    valid_candidates = []
    
    # Helper function: recursively determine maximum chain length along acyclic carbon atoms
    def get_max_chain_length(atom, visited):
        max_len = 0
        for nbr in atom.GetNeighbors():
            if nbr.GetIdx() in visited:
                continue
            if nbr.IsInRing() or nbr.GetAtomicNum() != 6:
                continue
            new_visited = visited | {nbr.GetIdx()}
            length = 1 + get_max_chain_length(nbr, new_visited)
            if length > max_len:
                max_len = length
        return max_len

    # Iterate over each ester candidate found by the SMARTS query.
    for match in ester_matches:
        # By the SMARTS "O[C](=O)", match[0] is the ester oxygen and match[1] is the carbonyl carbon.
        oxy_idx, carbonyl_idx = match[0], match[1]
        atom_oxy = mol.GetAtomWithIdx(oxy_idx)
        atom_caro = mol.GetAtomWithIdx(carbonyl_idx)
        
        # For the ester oxygen, find the substituent that is not the carbonyl carbon.
        oxy_neighbors = [nbr for nbr in atom_oxy.GetNeighbors() if nbr.GetIdx() != carbonyl_idx]
        if not oxy_neighbors:
            continue  # Ester oxygen without additional substituent is invalid.
        chol_candidate = oxy_neighbors[0]
        # The cholesterol part should be embedded in a ring.
        if not chol_candidate.IsInRing():
            continue
        
        # Check that the cholesterol candidate is part of a ring system in the nucleus.
        found_in_ring = False
        for ring in rings:
            if chol_candidate.GetIdx() in ring:
                found_in_ring = True
                break
        if not found_in_ring:
            continue
        
        # For the carbonyl carbon, identify the fatty acid substituent (neighbor not being the oxygen).
        fatty_neighbors = [nbr for nbr in atom_caro.GetNeighbors() if nbr.GetIdx() != oxy_idx]
        if not fatty_neighbors:
            continue
        fatty_candidate = fatty_neighbors[0]
        # The fatty acid chain should be acyclic.
        if fatty_candidate.IsInRing():
            continue
        
        # Determine the chain length by exploring connected acyclic carbon atoms.
        chain_length = get_max_chain_length(fatty_candidate, {fatty_candidate.GetIdx()})
        if chain_length < 5:
            continue
        
        # This candidate passes all tests.
        valid_candidates.append(match)
        
    if len(valid_candidates) == 0:
        return False, "No valid cholesteryl ester linkage found: ester group missing or not in proper cholesterol/fatty acid context"
    elif len(valid_candidates) > 1:
        return False, f"Multiple ({len(valid_candidates)}) ester linkages found; expected exactly one cholesteryl ester linkage"
    else:
        return True, "Molecule is classified as a cholesteryl ester with a steroid nucleus and a valid fatty acid chain"

# Example usage:
if __name__ == "__main__":
    # Test one example from the provided list: cholesteryl oleate
    test_smiles = "[H][C@@]1(CC[C@@]2([H])[C@]3([H])CC=C4C[C@H](CC[C@]4(C)[C@@]3([H])CC[C@]12C)OC(=O)CCCCCCC\\C=C/CCCCCCCC)"
    result, reason = is_cholesteryl_ester(test_smiles)
    print(result, reason)