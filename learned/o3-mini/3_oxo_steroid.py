"""
Classifies: CHEBI:47788 3-oxo steroid
"""
"""
Classifies: 3-oxo steroid
Definition: Any oxo steroid where an oxo substituent is located at position 3.
The function is_3_oxo_steroid takes a SMILES string as input and returns a tuple 
(bool, reason) indicating whether the molecule is a 3-oxo steroid and why.
"""

from rdkit import Chem

def is_3_oxo_steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo steroid based on its SMILES string.
    A 3-oxo steroid is defined as a steroid (typical tetracyclic nucleus: three 6-membered rings 
    and one 5-membered ring fused together) possessing a ketone group (C=O) on the A ring 
    (i.e. the six-membered ring that is not fused with the 5-membered ring).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-oxo steroid, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get all ring information (each ring is a tuple of atom indices).
    ring_info = mol.GetRingInfo().AtomRings()
    if not ring_info:
        return False, "No rings found in molecule, so not a steroid"

    # Convert rings to a list of sets of atom indices for easier manipulation.
    rings = [set(r) for r in ring_info]

    # Now group rings into fused ring systems.
    # Two rings are considered fused if they share at least one atom.
    n = len(rings)
    visited = [False] * n
    components = []
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

    # Look for a fused component that has exactly 4 rings with sizes 5,6,6,6.
    steroid_component = None
    for comp in components:
        sizes = sorted([len(rings[i]) for i in comp])
        if len(comp) == 4 and sizes == [5, 6, 6, 6]:
            steroid_component = comp
            break
    if steroid_component is None:
        return False, "Steroid nucleus not found (expected fused 4-ring system with rings of sizes 5,6,6,6)"

    # Identify the five-membered ring (assumed to be the D ring) and collect the six-membered rings.
    d_ring = None
    six_rings = []
    for i in steroid_component:
        ring_size = len(rings[i])
        if ring_size == 5:
            d_ring = rings[i]
        elif ring_size == 6:
            six_rings.append(rings[i])
    if d_ring is None or len(six_rings) != 3:
        return False, "Steroid nucleus ring pattern is not as expected"

    # Heuristically determine the A ring: it is the six-membered ring that does NOT share any 
    # atom with the five-membered (D ring).
    ringA = None
    for r in six_rings:
        if r.isdisjoint(d_ring):
            ringA = r
            break
    if ringA is None:
        return False, "Unable to identify the A ring of the steroid nucleus"

    # Now check for a ketone (C=O) on the identified A ring.
    # We iterate over the carbon atoms in the A ring and see if one has a double bond to oxygen.
    ketone_found = False
    for idx in ringA:
        atom = mol.GetAtomWithIdx(idx)
        if atom.GetAtomicNum() != 6:
            continue  # Only consider carbon atoms.
        for bond in atom.GetBonds():
            # Check if the bond is a double bond.
            if bond.GetBondType() == Chem.BondType.DOUBLE:
                # Identify the neighbor atom.
                neighbor = bond.GetOtherAtom(atom)
                if neighbor.GetAtomicNum() == 8:
                    ketone_found = True
                    break
        if ketone_found:
            break

    if not ketone_found:
        return False, "Steroid nucleus found but no ketone group detected on the A ring (position 3)"
    
    return True, "Steroid nucleus with a ketone group at position 3 identified"
    
# (Optional) If running as a script, one might include some tests, for example:
if __name__ == "__main__":
    test_smiles = [
        "[H][C@@]12CC(=O)CC[C@]1(C)[C@@]1([H])CC(=O)[C@]3(C)[C@]([H])(CC[C@@]3([H])[C@H](C)CCC(O)=O",  # 3,7,12-trioxo-5beta-cholanic acid
        "[H][C@@]12CC[C@]3([H])[C@]([H])(CC[C@]4(C)C=CC[C@@]34[H])[C@@]1(C)CCC(=O)C2"  # 5alpha-androst-16-en-3-one
    ]
    for s in test_smiles:
        result, reason = is_3_oxo_steroid(s)
        print(f"SMILES: {s}\nClassification: {result}\nReason: {reason}\n")