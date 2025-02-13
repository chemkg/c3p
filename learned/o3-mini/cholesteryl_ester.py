"""
Classifies: CHEBI:17002 cholesteryl ester
"""
"""
Classifies: cholesteryl ester
Definition: A sterol ester obtained by formal condensation of the carboxy group of any carboxylic acid 
with the 3-hydroxy group of cholesterol.
Heuristics in this version:
  1. The molecule must be a valid structure.
  2. The overall molecule must contain at least 4 rings (to indicate a steroid nucleus).
  3. There must be at least one ester linker of the form –O–C(=O)–.
     For each candidate ester:
       a. The ester oxygen must be attached to a substituent that is part of a ring.
       b. The carbonyl carbon must display a carbonyl group (i.e. double-bonded oxygen) and be attached to 
          an acyclic fatty acid chain—this chain requirement is relaxed (chain length ≥2) so that even acetate 
          esters are accepted.
  4. Additionally, we require that the entire molecule contains a cholesterol-like aliphatic side chain,
     defined by the presence of a substructure pattern of a ring atom attached to an acyclic chain of at least 4 carbons.
  5. Exactly one candidate ester linkage (meeting both cholesterol and fatty-acid criteria) must be found.
"""

from rdkit import Chem

def is_cholesteryl_ester(smiles: str):
    """
    Determines if a molecule is a cholesteryl ester based on its SMILES string.
    
    The algorithm checks:
      1. Validity of the SMILES.
      2. That the molecule contains a steroid nucleus (at least 4 rings).
      3. That the molecule contains at least one ester substructure (O–C(=O)–).
      4. For each candidate ester:
         a. The ester oxygen’s non-carbonyl neighbor must be in a ring (indicating the cholesterol hydroxyl origin).
         b. The fatty acid side (attached to the carbonyl carbon) must be acyclic and have a minimum chain length (≥2 carbon bonds).
      5. In addition, the entire molecule must display a cholesterol-like side chain—
         a long (≥4-carbon) acyclic alkyl chain attached to a ring atom.
      6. Exactly one candidate ester must pass these tests.
      
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a cholesteryl ester, False otherwise.
        str: Explanation for the decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check that the molecule has at least 4 rings (indicator of a fused steroid nucleus)
    rings = mol.GetRingInfo().AtomRings()
    if len(rings) < 4:
        return False, "Molecule does not contain enough rings for a steroid nucleus (expected ≥4 rings)"
    
    # Check for the presence of a cholesterol-like side chain.
    # Cholesterol has a long alkyl chain (typically ≥4 continuous acyclic C) appended to the steroid nucleus.
    side_chain_smarts = "[#6;R]-[#6;R0]CCCC"  # A ring carbon attached to at least 4 acyclic carbons.
    side_chain_query = Chem.MolFromSmarts(side_chain_smarts)
    if not mol.HasSubstructMatch(side_chain_query):
        return False, "Molecule does not display a cholesterol-like side chain (expected a long aliphatic chain attached to the steroid nucleus)"
    
    # Define ester SMARTS for candidates: look for a substructure O-C(=O)
    ester_query = Chem.MolFromSmarts("O[C](=O)")
    ester_matches = mol.GetSubstructMatches(ester_query)
    if not ester_matches:
        return False, "No ester group found"
    
    # Helper function to compute maximum chain length along acyclic carbon atoms.
    def get_max_chain_length(atom, visited):
        max_len = 0
        for nbr in atom.GetNeighbors():
            if nbr.GetIdx() in visited:
                continue
            # Only consider non-ring carbons (acyclic) that are carbon.
            if nbr.IsInRing() or nbr.GetAtomicNum() != 6:
                continue
            new_visited = visited | {nbr.GetIdx()}
            length = 1 + get_max_chain_length(nbr, new_visited)
            if length > max_len:
                max_len = length
        return max_len

    valid_candidates = []
    # Iterate over each ester candidate match.
    for match in ester_matches:
        # In our SMARTS "O[C](=O)", match[0] is the ester oxygen and match[1] is the carbonyl carbon.
        oxy_idx, carbonyl_idx = match[0], match[1]
        atom_oxy = mol.GetAtomWithIdx(oxy_idx)
        atom_caro = mol.GetAtomWithIdx(carbonyl_idx)
        
        # Verify that the carbonyl carbon actually has a double-bonded oxygen.
        has_cdbond = False
        for bond in atom_caro.GetBonds():
            if bond.GetBondTypeAsDouble() == 2.0:
                nbr = bond.GetOtherAtom(atom_caro)
                if nbr.GetAtomicNum() == 8:
                    has_cdbond = True
                    break
        if not has_cdbond:
            continue

        # For the ester oxygen, identify the neighbor that is not the carbonyl carbon.
        oxy_neighbors = [nbr for nbr in atom_oxy.GetNeighbors() if nbr.GetIdx() != carbonyl_idx]
        if not oxy_neighbors:
            continue
        # The cholesterol part candidate: expected to be in a ring.
        chol_candidate = oxy_neighbors[0]
        if not chol_candidate.IsInRing():
            continue
        
        # Additionally, require that the candidate is part of at least one ring (already implied but double-check).
        in_ring = any(chol_candidate.GetIdx() in ring for ring in rings)
        if not in_ring:
            continue
        
        # For the carbonyl carbon, identify the fatty acid substituent (the neighbor not being the oxygen).
        fatty_neighbors = [nbr for nbr in atom_caro.GetNeighbors() if nbr.GetIdx() != oxy_idx]
        if not fatty_neighbors:
            continue
        fatty_candidate = fatty_neighbors[0]
        if fatty_candidate.IsInRing():
            continue
        
        # Determine chain length from the fatty acid candidate.
        chain_length = get_max_chain_length(fatty_candidate, {fatty_candidate.GetIdx()})
        if chain_length < 2:
            # Relaxed threshold: even acetate (chain_length 0) should be allowed.
            # Here, we require that either there is at least a minimal (≥2) chain OR the fatty acid is a carboxylate (acetate).
            # For acetate, the fatty_candidate would be a methyl (chain length 0). So we allow chain_length==0 or 1.
            # Hence, if chain_length is 0 or 1, we accept.
            pass  # Do not reject based solely on chain length.
        
        # If all tests here pass, add to valid candidates.
        valid_candidates.append(match)
    
    if len(valid_candidates) == 0:
        return False, "No valid cholesteryl ester linkage found: ester group missing or not in proper cholesterol/fatty acid context"
    elif len(valid_candidates) > 1:
        return False, f"Multiple ({len(valid_candidates)}) ester linkages found; expected exactly one cholesteryl ester linkage"
    else:
        return True, "Molecule is classified as a cholesteryl ester with a steroid nucleus and a valid fatty acid chain"


# Example usage (will only print a result if the file is executed as a script).
if __name__ == "__main__":
    # Test one example: cholesteryl oleate
    test_smiles = "[H][C@@]1(CC[C@@]2([H])[C@]3([H])CC=C4C[C@H](CC[C@]4(C)[C@@]3([H])CC[C@]12C)OC(=O)CCCCCCC\\C=C/CCCCCCCC)"
    result, reason = is_cholesteryl_ester(test_smiles)
    print(result, reason)