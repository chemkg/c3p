"""
Classifies: CHEBI:15693 aldose
"""
"""
Classifies: Aldose – Aldehydic parent sugars (polyhydroxy aldehydes H[CH(OH)]nC(=O)H, n >= 2)
and their intramolecular hemiacetals.
Examples include beta-D-idopyranose, D-ribofuranose, L-erythrose, beta-D-glucose, etc.
"""

from rdkit import Chem

def is_aldose(smiles: str):
    """
    Determines if a molecule is an aldose based on its SMILES string.
    Aldoses can either be open-chain aldehydes with a fully hydroxylated backbone or
    cyclic hemiacetals having a 5- or 6-membered ring containing exactly one oxygen and
    appropriate exocyclic hydroxyl (–OH or CH2OH) decorations on each ring carbon.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule appears to be an aldose, False otherwise.
        str: Detailed reason for the classification decision.
    """
    # Parse the SMILES string to an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count carbon atoms in the molecule.
    carbons = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    nC = len(carbons)
    if nC < 3 or nC > 7:
        return False, f"Number of carbons is {nC}, which is outside the typical range for aldoses (3-7)"

    # Helper function: determines if an oxygen behaves as a hydroxyl (-OH).
    def is_hydroxyl(oxygen):
        if oxygen.GetAtomicNum() != 8:
            return False
        # A hydroxyl oxygen should have at least one hydrogen neighbor.
        return any(nb.GetAtomicNum() == 1 for nb in oxygen.GetNeighbors())

    # Helper for cyclic sugars:
    # Check if a given ring carbon has an exocyclic hydroxyl,
    # either directly or via a CH2OH group.
    def has_exocyclic_hydroxyl(atom, ring_indices):
        # Check each neighbor that is not part of the ring.
        for nb in atom.GetNeighbors():
            if nb.GetIdx() in ring_indices:
                continue
            # Direct attachment to a hydroxyl oxygen.
            if nb.GetAtomicNum() == 8 and is_hydroxyl(nb):
                return True
            # Alternatively, if a neighboring carbon (e.g., CH2) eventually bears an -OH.
            if nb.GetAtomicNum() == 6:
                for nb2 in nb.GetNeighbors():
                    if nb2.GetIdx() == atom.GetIdx() or nb2.GetIdx() in ring_indices:
                        continue
                    if nb2.GetAtomicNum() == 8 and is_hydroxyl(nb2):
                        return True
        return False

    # Attempt 1: Open-chain aldose detection (only if the molecule is acyclic).
    if mol.GetRingInfo().NumRings() == 0:
        # Use SMARTS to detect a terminal aldehyde: a carbonyl (C=O) carbon with one hydrogen.
        aldehyde_smarts = Chem.MolFromSmarts("[CX3H1](=O)")
        aldehyde_matches = mol.GetSubstructMatches(aldehyde_smarts)
        for match in aldehyde_matches:
            aldehyde_atom = mol.GetAtomWithIdx(match[0])
            # Check that the aldehyde carbon is terminal (attached to exactly one other carbon).
            carbon_neighbors = [nb for nb in aldehyde_atom.GetNeighbors() if nb.GetAtomicNum() == 6]
            if len(carbon_neighbors) != 1:
                continue  # Not a terminal aldehyde.
            # Traverse the carbon chain from the aldehyde carbon.
            chain = set()
            def dfs(atom):
                if atom.GetIdx() in chain:
                    return
                chain.add(atom.GetIdx())
                for nb in atom.GetNeighbors():
                    if nb.GetAtomicNum() == 6 and nb.GetIdx() not in chain:
                        dfs(nb)
            dfs(aldehyde_atom)
            if len(chain) != nC: 
                continue  # The carbon chain is not contiguous.
            # For every carbon (except the terminal aldehyde), verify there is an exocyclic -OH.
            open_chain_valid = True
            for cid in chain:
                if cid == aldehyde_atom.GetIdx():
                    continue
                carbon_atom = mol.GetAtomWithIdx(cid)
                found_hydroxyl = False
                for nb in carbon_atom.GetNeighbors():
                    if nb.GetIdx() in chain:
                        continue
                    if nb.GetAtomicNum() == 8 and is_hydroxyl(nb):
                        found_hydroxyl = True
                        break
                    if nb.GetAtomicNum() == 6:
                        # Check if this substituent carbon carries an -OH.
                        for nb2 in nb.GetNeighbors():
                            if nb2.GetIdx() == carbon_atom.GetIdx() or nb2.GetIdx() in chain:
                                continue
                            if nb2.GetAtomicNum() == 8 and is_hydroxyl(nb2):
                                found_hydroxyl = True
                                break
                        if found_hydroxyl:
                            break
                if not found_hydroxyl:
                    open_chain_valid = False
                    break
            if open_chain_valid:
                return True, ("Open-chain aldose pattern detected with terminal aldehyde "
                              "and fully hydroxylated carbon backbone")
        return False, "No open-chain aldose pattern detected"

    # Attempt 2: Cyclic (hemiacetal) aldose detection.
    rings = mol.GetRingInfo().AtomRings()
    for ring in rings:
        if len(ring) not in (5, 6):
            continue  # Only consider 5- or 6-membered rings.
        ring_indices = set(ring)
        # For an aldose ring, expect exactly one oxygen (the ring oxygen).
        oxygen_in_ring = [idx for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8]
        if len(oxygen_in_ring) != 1:
            continue
        # Check that the ring (excluding the ring oxygen) is composed of carbons.
        ring_carbons = [idx for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6]
        if len(ring_carbons) != len(ring) - 1:
            continue
        # Now ensure that each ring carbon has an exocyclic hydroxyl decoration.
        valid_ring = True
        for idx in ring_carbons:
            atom = mol.GetAtomWithIdx(idx)
            if not has_exocyclic_hydroxyl(atom, ring_indices):
                valid_ring = False
                break
        if valid_ring:
            return True, (f"Cyclic hemiacetal sugar pattern detected in a {len(ring)}-membered "
                          "ring with appropriate exocyclic hydroxyl decorations")
    
    # If neither open-chain nor cyclic patterns were detected, classify as non-aldose.
    return False, ("No aldose pattern detected: neither open-chain terminal aldehyde with a "
                   "fully hydroxylated backbone nor a cyclic hemiacetal sugar pattern was found")

# Example test cases (if running as a main file).
if __name__ == "__main__":
    test_cases = [
        # Open-chain aldoses:
        ("[H]C(=O)[C@@H](O)[C@H](O)CO", "aldehydo-L-ribose (open-chain)"),
        ("[H]C(=O)[C@H](O)CO", "D-glyceraldehyde (open-chain)"),
        ("[H][C@](O)(CO)[C@]([H])(O)C=O", "L-erythrose (open-chain)"),
        # Cyclic aldoses:
        ("O1[C@@H]([C@H](O)[C@@H](O)[C@H](O)[C@@H]1O)CO", "beta-D-idopyranose (cyclic)"),
        ("O1[C@@H](O)[C@@H](O)[C@@H](O)[C@H](O)C1", "beta-D-lyxopyranose (cyclic)"),
        ("OC[C@H]1OC(O)[C@H](O)[C@@H]1O", "D-ribofuranose (cyclic)"),
        ("O1[C@H]([C@H](O)[C@H](O)C1O)[C@H](O)CO", "D-talofuranose (cyclic)"),
        ("OC[C@@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O", "beta-D-glucose (cyclic)")
    ]
    for smi, name in test_cases:
        result, reason = is_aldose(smi)
        print(f"SMILES: {smi} | {name} -> {result} ; Reason: {reason}")