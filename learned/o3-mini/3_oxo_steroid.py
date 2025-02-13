"""
Classifies: CHEBI:47788 3-oxo steroid
"""
"""
Classifies: 3-oxo steroid
Definition: Any oxo steroid where an oxo substituent is located at position 3.
This algorithm detects a steroid nucleus by trying to find a fused tetracyclic system (4 rings)
with relaxed ring size constraints (one small [≤5] and other rings between 5 and 7 atoms) and a total of
roughly 16–24 carbon atoms in the fused core. Then, among the six-membered rings in that core, it looks for a
ketone group (C=O) with a terminal oxygen, which is assumed to be at the 3‐position.
"""

from rdkit import Chem

def is_3_oxo_steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo steroid based on its SMILES.
    A 3-oxo steroid is defined as a steroid having a fused tetracyclic nucleus (typically a cyclopentanoperhydrophenanthrene
    system) and possessing a ketone group in one of its six-membered rings (a proxy for a 3-oxo substituent).
    
    The algorithm:
      1. Parses the molecule and obtains all rings.
      2. Groups rings into fused clusters.
      3. Selects candidate clusters that have four fused rings and that meet relaxed size criteria:
           - At least one ring must be “small” (≤5 members) and the others between 5 and 7 atoms.
           - The union of all atoms in the candidate has a carbon count between 16 and 24.
      4. Within the candidate, it looks for at least one six-membered ring containing a terminal ketone.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        (bool, str): True with explanation if molecule is classified as 3-oxo steroid, or False with reason.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get ring information
    ring_info = mol.GetRingInfo().AtomRings()
    if not ring_info:
        return False, "No rings found in molecule, so not a steroid"

    # Convert each ring (tuple of indices) to a set for easier processing.
    rings = [set(r) for r in ring_info]

    # Group rings into fused clusters. Rings are fused if they share at least one atom.
    n = len(rings)
    visited = [False] * n
    fused_clusters = []  # each cluster will be a set of ring indices
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
                    # if rings share any atom then they are fused
                    if rings[current].intersection(rings[j]):
                        comp.add(j)
                        stack.append(j)
                        visited[j] = True
        fused_clusters.append(comp)

    candidate_cluster = None
    # We expect steroid nucleus as a fused system of four rings.
    for comp in fused_clusters:
        if len(comp) != 4:
            continue
        # get the ring sizes in this cluster:
        sizes = [len(rings[i]) for i in comp]
        # relaxed requirement: one ring is small (<=5) and others have size between 5 and 7.
        if not any(s <= 5 for s in sizes):
            continue
        if not all(5 <= s <= 7 for s in sizes):
            continue

        # compute the union of all atoms in the candidate core:
        core_atoms = set()
        for i in comp:
            core_atoms = core_atoms.union(rings[i])
        # count carbon atoms in the fused core
        core_carbon_count = sum(1 for idx in core_atoms if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
        if core_carbon_count < 16 or core_carbon_count > 24:
            continue

        candidate_cluster = comp
        break

    if candidate_cluster is None:
        return False, "Steroid nucleus not found (expected fused tetracyclic system with relaxed ring size and carbon count criteria)"

    # Get union of atoms in candidate cluster:
    core_atoms = set()
    for i in candidate_cluster:
        core_atoms = core_atoms.union(rings[i])
    
    # Now, among the rings in the candidate, collect those that are six-membered.
    six_membered_rings = [rings[i] for i in candidate_cluster if len(rings[i]) == 6]
    if not six_membered_rings:
        return False, "Fused ring system found but no six-membered ring present for ketone detection"

    # Look for ketone in one of the six-membered rings.
    # A ketone is defined as a carbon atom in the ring that is double bonded to an oxygen (atomic num 8)
    # where the oxygen is terminal (only one non-hydrogen neighbor).
    ketone_found = False
    ketone_detail = ""
    for ring in six_membered_rings:
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() != 6:
                continue
            for bond in atom.GetBonds():
                # we check that the bond is a double bond.
                if bond.GetBondType() == Chem.BondType.DOUBLE:
                    neighbor = bond.GetOtherAtom(atom)
                    if neighbor.GetAtomicNum() == 8:
                        # Check that the oxygen is terminal (only one non-hydrogen neighbor)
                        nonH_neighbors = [nbr for nbr in neighbor.GetNeighbors() if nbr.GetAtomicNum() != 1]
                        if len(nonH_neighbors) == 1:
                            ketone_found = True
                            ketone_detail = f"Ketone (C=O) found on atom index {idx} in one six-membered ring"
                            break
            if ketone_found:
                break
        if ketone_found:
            break

    if not ketone_found:
        return False, "Steroid nucleus found but no terminal ketone group in any six-membered ring (expected at position 3)"

    return True, f"Steroid nucleus with {core_carbon_count} core carbons and a ketone group detected: {ketone_detail}"


# Optional main block for testing examples.
if __name__ == "__main__":
    test_smiles = [
        # True positives (should be classified as 3-oxo steroids):
        "[H][C@@]12CC(=O)CC[C@]1(C)[C@@]1([H])CC(=O)[C@]3(C)[C@]([H])(CC[C@@]3([H])[C@H](C)CCC(O)=O",  # 3,7,12-trioxo-5beta-cholanic acid
        "[H][C@@]12CC[C@]3([H])[C@]([H])(CC[C@]4(C)C=CC[C@@]34[H])[C@@]1(C)CCC(=O)C2",  # 5alpha-androst-16-en-3-one
        # False positives and negatives can be further tested with the provided SMILES in the prompt.
    ]
    
    for s in test_smiles:
        result, reason = is_3_oxo_steroid(s)
        print("SMILES:", s)
        print("Classification:", result)
        print("Reason:", reason)
        print("-----")