"""
Classifies: CHEBI:68472 pyrimidine deoxyribonucleoside
"""
"""
Classifies: A deoxyribonucleoside containing a pyrimidine base.
A pyrimidine deoxyribonucleoside contains a five‐membered deoxyribose sugar (a furanose ring missing the 2'-OH) 
that is linked via a glycosidic bond to a pyrimidine base.
"""
from rdkit import Chem

def is_pyrimidine_deoxyribonucleoside(smiles: str):
    """
    Determines if a molecule is a pyrimidine deoxyribonucleoside.
    The method uses a heuristic based on:
      - detecting a pyrimidine base via a generic SMARTS "c1cncnc1"
      - detecting a five-membered sugar ring (furanose) that contains exactly one ring oxygen and 
        exactly 2 exocyclic oxygen atoms attached to ring carbon atoms [indicative of a deoxyribose sugar].
      - confirming that one atom from the sugar ring is directly connected to an atom in the pyrimidine base.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is classified as a pyrimidine deoxyribonucleoside, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Step 1: Look for a pyrimidine base using a generic SMARTS pattern.
    # The pattern "c1cncnc1" matches a six-membered aromatic ring alternating between C and N.
    pyrimidine_smarts = "c1cncnc1"
    pyrimidine_query = Chem.MolFromSmarts(pyrimidine_smarts)
    if not mol.HasSubstructMatch(pyrimidine_query):
        return False, "Pyrimidine base not found"
    pyrimidine_matches = mol.GetSubstructMatches(pyrimidine_query)
    # Take the first match as a candidate pyrimidine (set of atom indices).
    pyr_atom_indices = set(pyrimidine_matches[0])

    # Step 2: Look for a deoxyribose sugar moiety.
    # We need a five-membered ring (furanose) that contains exactly one oxygen atom (in-ring)
    # and where exocyclic oxygen atoms (attached to ring carbon atoms and carrying at least one hydrogen)
    # number exactly 2 (loss of the 2'-OH gives deoxyribose, whereas ribose would have 3).
    ring_info = mol.GetRingInfo()
    candidate_sugar = None
    for ring in ring_info.AtomRings():
        if len(ring) == 5:
            # Check how many atoms in the ring are oxygens.
            n_ring_oxy = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetSymbol() == "O")
            if n_ring_oxy != 1:
                continue  # Not a furanose ring
            # Count exocyclic oxygen attachments on ring carbons (exclude the ring oxygen already counted)
            exo_oxy_count = 0
            for idx in ring:
                atom = mol.GetAtomWithIdx(idx)
                # Only consider carbon atoms (since the ring oxygen is part of the ring)
                if atom.GetSymbol() != "C":
                    continue
                for nbr in atom.GetNeighbors():
                    if nbr.GetIdx() not in ring and nbr.GetSymbol() == "O":
                        # Check if this oxygen has at least one hydrogen attached (i.e. is -OH or part of CH2OH)
                        if any(neigh.GetAtomicNum() == 1 for neigh in nbr.GetNeighbors()):
                            exo_oxy_count += 1
            # For a deoxyribose in a nucleoside, we expect exactly 2 exocyclic oxygen attachments:
            # one from the 3'-OH and one from the 5'-CH2OH (the 1' position links to the base).
            if exo_oxy_count != 2:
                continue

            # Step 3: Verify that this sugar ring is linked to the pyrimidine base by a glycosidic bond.
            sugar_linked = False
            for idx in ring:
                atom = mol.GetAtomWithIdx(idx)
                for nbr in atom.GetNeighbors():
                    if nbr.GetIdx() in pyr_atom_indices:
                        sugar_linked = True
                        break
                if sugar_linked:
                    break

            if not sugar_linked:
                continue

            # If we reach here, we found a candidate sugar ring.
            candidate_sugar = ring
            break

    if candidate_sugar is None:
        return False, "Deoxyribose sugar moiety not found"

    return True, "Molecule contains a pyrimidine base linked to a deoxyribose sugar"


# For testing purposes, you can try calling the function with one of the provided SMILES.
if __name__ == "__main__":
    # Example: 2'-deoxyuridine
    test_smiles = "OC[C@H]1O[C@H](C[C@@H]1O)n1ccc(=O)[nH]c1=O"
    result, reason = is_pyrimidine_deoxyribonucleoside(test_smiles)
    print(result, reason)