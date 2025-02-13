"""
Classifies: CHEBI:47788 3-oxo steroid
"""
"""
Classifies: 3-oxo steroid
Definition: Any oxo steroid where an oxo substituent is located at position 3.
This algorithm first finds fused ring systems in the molecule and then looks for a
component that has four fused rings with sizes 5,6,6,6 and a core carbon count between 17 and 21.
Then, among the six membered rings in that fused core, it looks for a ketone group (C=O).
The function is_3_oxo_steroid takes a SMILES string as input and returns a tuple (bool, reason)
indicating whether the molecule is a 3-oxo steroid and why.
"""

from rdkit import Chem

def is_3_oxo_steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo steroid based on its SMILES string.
    A 3-oxo steroid is defined as a steroid (with a fused tetracyclic nucleus: one five-membered ring
    and three six-membered rings that together contain between roughly 17 and 21 carbons) possessing a ketone group (C=O)
    on one of the six-membered rings (i.e. in the A ring region).
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule is a 3-oxo steroid, False otherwise
        str: Explanation for the classification decision
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get ring information from the molecule.
    ring_info = mol.GetRingInfo().AtomRings()
    if not ring_info:
        return False, "No rings found in molecule, so not a steroid"

    # Convert each ring (tuple of atom indices) to a set for easier processing.
    rings = [set(r) for r in ring_info]

    # Group rings into fused components: rings are fused if they share at least one atom.
    n = len(rings)
    visited = [False] * n
    components = []  # each component is a set of ring indices.
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
                    if rings[current].intersection(rings[j]):
                        comp.add(j)
                        stack.append(j)
                        visited[j] = True
        components.append(comp)

    candidate_component = None
    # Look for a component with exactly 4 rings where the ring sizes (sorted) are [5,6,6,6].
    for comp in components:
        ring_sizes = sorted([len(rings[i]) for i in comp])
        if len(comp) == 4 and ring_sizes == [5, 6, 6, 6]:
            candidate_component = comp
            break
    if candidate_component is None:
        return False, "Steroid nucleus not found (expected four fused rings with sizes 5,6,6,6)"

    # Get the union of all atom indices in the candidate component.
    core_atom_indices = set()
    for i in candidate_component:
        core_atom_indices = core_atom_indices.union(rings[i])
    
    # Count carbon atoms in the steroid core.
    core_carbon_count = sum(1 for idx in core_atom_indices if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
    if core_carbon_count < 17 or core_carbon_count > 21:
        return False, f"Fused ring system found but with {core_carbon_count} core carbons (expected between 17 and 21)"

    # Among the rings in the candidate component, collect the six-membered rings.
    six_membered_rings = [rings[i] for i in candidate_component if len(rings[i]) == 6]
    if len(six_membered_rings) != 3:
        return False, "Fused ring system present but the number of six-membered rings is not equal to 3"

    # Now, for each six-membered ring in the core, check for a ketone.
    # A ketone group is defined as a carbon (atomic number 6) in the ring that has
    # at least one double bond to an oxygen atom (atomic number 8) and that oxygen should be terminal.
    ketone_found = False
    ketone_details = ""
    for ring in six_membered_rings:
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() != 6:
                continue
            for bond in atom.GetBonds():
                # Check that the bond is a double bond.
                if bond.GetBondType() == Chem.BondType.DOUBLE:
                    neighbor = bond.GetOtherAtom(atom)
                    if neighbor.GetAtomicNum() == 8:
                        # Check that the oxygen is terminal (only one non-hydrogen neighbor)
                        nonH_neighbors = [nbr for nbr in neighbor.GetNeighbors() if nbr.GetAtomicNum() != 1]
                        if len(nonH_neighbors) == 1:
                            ketone_found = True
                            ketone_details = f"Ketone (C=O) found on atom index {idx} in a six-membered ring"
                            break
            if ketone_found:
                break
        if ketone_found:
            break

    if not ketone_found:
        return False, "Steroid nucleus core found but no ketone group detected on any six-membered ring (expected on position 3)"
    
    return True, f"Steroid nucleus with {core_carbon_count} core carbons and a ketone group detected: {ketone_details}"


# Optional main block for testing examples.
if __name__ == "__main__":
    test_smiles = [
        # Examples that should be classified as 3-oxo steroids:
        "[H][C@@]12CC(=O)CC[C@]1(C)[C@@]1([H])CC(=O)[C@]3(C)[C@]([H])(CC[C@@]3([H])[C@H](C)CCC(O)=O",  # 3,7,12-trioxo-5beta-cholanic acid
        "[H][C@@]12CC[C@]3([H])[C@]([H])(CC[C@]4(C)C=CC[C@@]34[H])[C@@]1(C)CCC(=O)C2"  # 5alpha-androst-16-en-3-one
    ]
    for s in test_smiles:
        result, reason = is_3_oxo_steroid(s)
        print("SMILES:", s)
        print("Classification:", result)
        print("Reason:", reason)
        print("-----")