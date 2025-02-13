"""
Classifies: CHEBI:47788 3-oxo steroid
"""
"""
Classifies: 3-oxo steroid
Definition: Any oxo steroid where an oxo substituent is located at position 3.
This version attempts to detect a steroid nucleus by finding a fused cluster
of exactly four rings (ideally three six‐membered and one five‐membered), where the union
of the rings gives between 16 and 18 carbon atoms (a proxy for a 17-carbon nucleus).
Then the algorithm restricts ketone detection to carbons within the nucleus that belong 
to a six-membered ring and are double-bonded to a terminal oxygen.
"""

from rdkit import Chem

def is_3_oxo_steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo steroid based on its SMILES string.
    The algorithm:
      1. Parses the molecule.
      2. Retrieves all rings and groups those that are fused (sharing one or more atoms).
      3. For each fused cluster of exactly 4 rings the union of ring atoms is considered.
         If the cluster contains between 16 and 18 carbon atoms (a proxy for the steroid core)
         and at least one of the rings is six-membered, it is retained.
      4. Within any six-membered ring of the candidate fused cluster, the function searches
         for a carbon that bears a double bond to oxygen. The oxygen must be terminal
         (i.e. it only bonds to that carbon aside from any hydrogen atoms). This is taken as
         evidence of a 3-oxo substituent.
         
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        (bool, str): True with explanation if molecule is classified as a 3-oxo steroid,
                     False along with an explanation otherwise.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Retrieve ring information from the molecule.
    ring_info = mol.GetRingInfo().AtomRings()
    if not ring_info:
        return False, "No rings found in molecule, so not a steroid"

    # Convert each ring from tuple to set for easier processing.
    rings = [set(r) for r in ring_info]
    n = len(rings)
    visited = [False] * n
    fused_clusters = []  # each cluster is a set of ring indices

    # Group rings into fused clusters (rings sharing at least one atom are fused).
    for i in range(n):
        if visited[i]:
            continue
        comp = set([i])
        stack = [i]
        while stack:
            current = stack.pop()
            visited[current] = True
            for j in range(n):
                if not visited[j]:
                    if rings[current].intersection(rings[j]):  # fused if share an atom
                        comp.add(j)
                        stack.append(j)
                        visited[j] = True
        fused_clusters.append(comp)

    candidate_cluster = None
    candidate_core_atoms = None

    # Look through fused clusters for one that is exactly 4 rings (typical steroid nucleus).
    for comp in fused_clusters:
        if len(comp) != 4:
            continue
        # Calculate union of all atoms in the cluster.
        core_atoms = set()
        for idx in comp:
            core_atoms = core_atoms.union(rings[idx])
        # Count how many carbons are in the core.
        core_carbon_count = sum(1 for idx in core_atoms if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
        # We expect a steroid nucleus to have about 17 carbons (+/- 1). Adjust as needed.
        if core_carbon_count < 16 or core_carbon_count > 18:
            continue

        # Require that at least one ring in the candidate is six-membered.
        six_membered_found = any(len(rings[i]) == 6 for i in comp)
        if not six_membered_found:
            continue

        candidate_cluster = comp
        candidate_core_atoms = core_atoms
        break

    if candidate_cluster is None:
        return False, "Steroid nucleus not found (expected fused tetracyclic system with 16-18 core carbons)"

    # Among the candidate cluster, extract rings that are six-membered.
    six_membered_rings = [rings[i] for i in candidate_cluster if len(rings[i]) == 6]
    if not six_membered_rings:
        return False, "Fused ring system found but no six-membered ring present for ketone detection"

    # Look for a ketone (C=O) on an atom that belongs to the nucleus.
    ketone_found = False
    ketone_detail = ""
    # Iterate over each six-membered ring in the candidate fused cluster.
    for ring in six_membered_rings:
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            # Focus only on carbon atoms that are part of the candidate nucleus.
            if atom.GetAtomicNum() != 6:
                continue
            # Check all bonds from the carbon.
            for bond in atom.GetBonds():
                # Only consider double bonds.
                if bond.GetBondType() == Chem.BondType.DOUBLE:
                    neighbor = bond.GetOtherAtom(atom)
                    if neighbor.GetAtomicNum() == 8:
                        # Verify that the oxygen is terminal, i.e. it has only one non-hydrogen neighbor.
                        nonH_neighbors = [nbr for nbr in neighbor.GetNeighbors() if nbr.GetAtomicNum() != 1]
                        if len(nonH_neighbors) == 1:
                            ketone_found = True
                            ketone_detail = f"Ketone (C=O) found on atom index {idx} in a six-membered ring"
                            break
            if ketone_found:
                break
        if ketone_found:
            break

    if not ketone_found:
        return False, "Steroid nucleus found but no terminal ketone group detected in any six-membered ring (expected at position 3)"

    return True, f"Steroid nucleus with {sum(1 for idx in candidate_core_atoms if mol.GetAtomWithIdx(idx).GetAtomicNum()==6)} core carbons and a ketone group detected: {ketone_detail}"


# Optional main block for testing.
if __name__ == "__main__":
    test_cases = [
        # True positive examples:
        # 3,7,12-trioxo-5beta-cholanic acid
        "[H][C@@]12CC(=O)CC[C@]1(C)[C@@]1([H])CC(=O)[C@]3(C)[C@]([H])(CC[C@@]3([H])[C@H](C)CCC(O)=O",
        # 5alpha-androst-16-en-3-one
        "[H][C@@]12CC[C@]3([H])[C@]([H])(CC[C@]4(C)C=CC[C@@]34[H])[C@@]1(C)CCC(=O)C2",
        # False positive example (cucurbitacin D should not classify):
        "[H][C@@]12C[C@H](O)C(=O)C(C)(C)C1=CC[C@@]1([H])[C@]3(C)C[C@@H](O)[C@]([H])([C@@](C)(O)C(=O)\\C=C\\C(C)(C)O)[C@@]3(C)CC(=O)[C@@]21C"
    ]
    
    for s in test_cases:
        result, reason = is_3_oxo_steroid(s)
        print("SMILES:", s)
        print("Classification:", result)
        print("Reason:", reason)
        print("-----")