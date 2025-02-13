"""
Classifies: CHEBI:71971 neoflavonoid
"""
#!/usr/bin/env python
"""
Classifies: neoflavonoid
Definition: A neoflavonoid is any 1-benzopyran with an aryl substituent at position 4.
Heuristic:
  1. Identify candidate pyran rings: six-membered rings containing exactly one oxygen.
  2. Identify candidate benzene rings: six-membered rings where every atom is aromatic.
  3. Choose a pyran ring that is fused with a benzene ring (sharing at least two atoms).
  4. In that fused pyran ring, locate the unique oxygen and compute the bond distances along the ring.
     The candidate "position 4" is the atom at bond distance 3 from the oxygen, provided it is not part of the fused benzene.
  5. Verify that the candidate position bears exactly one substituent (via a single bond) that is part of an independent aromatic ring.
"""

from rdkit import Chem

def is_neoflavonoid(smiles: str):
    """
    Determines if a molecule is a neoflavonoid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a neoflavonoid, False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    try:
        Chem.SanitizeMol(mol)
    except Exception as e:
        return False, "Sanitization error: " + str(e)
    
    # Retrieve ring information
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    
    # Step 1: Identify candidate pyran rings (six-membered rings with exactly one oxygen atom)
    pyran_rings = []
    for ring in atom_rings:
        if len(ring) == 6:
            oxy_count = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8)
            if oxy_count == 1:
                pyran_rings.append(set(ring))
    if not pyran_rings:
        return False, "No six-membered oxygen-containing ring (pyran) found."
    
    # Step 2: Identify candidate benzene rings: aromatic six-membered rings.
    benzene_rings = []
    for ring in atom_rings:
        if len(ring) == 6:
            # Ensure every atom in the ring is aromatic
            if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
                benzene_rings.append(set(ring))
    if not benzene_rings:
        return False, "No aromatic benzene ring found, required for benzopyran core."
    
    # Step 3: Look for a pyran ring that is fused to a benzene ring (sharing at least two atoms)
    fused_pyran = None
    fused_benzene = None
    for pyran in pyran_rings:
        for benzene in benzene_rings:
            shared = pyran.intersection(benzene)
            if len(shared) >= 2:
                fused_pyran = pyran
                fused_benzene = benzene
                break
        if fused_pyran is not None:
            break
    if fused_pyran is None:
        return False, "No fused benzopyran core (pyran fused to benzene) detected."
    
    # Step 4: In the fused pyran, locate the unique oxygen atom.
    oxygen_atoms = [idx for idx in fused_pyran if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8]
    if len(oxygen_atoms) != 1:
        return False, "Fused pyran ring does not have exactly one oxygen atom."
    oxy_idx = oxygen_atoms[0]
    
    # Build connectivity graph within the fused pyran ring
    connectivity = { idx: set() for idx in fused_pyran }
    for idx in fused_pyran:
        atom = mol.GetAtomWithIdx(idx)
        for nbr in atom.GetNeighbors():
            nbr_idx = nbr.GetIdx()
            if nbr_idx in fused_pyran:
                connectivity[idx].add(nbr_idx)
    
    # Use BFS to determine the distance (number of bonds) from the oxygen atom within the ring.
    distances = { oxy_idx: 0 }
    queue = [oxy_idx]
    while queue:
        current = queue.pop(0)
        for nbr in connectivity[current]:
            if nbr not in distances:
                distances[nbr] = distances[current] + 1
                queue.append(nbr)
    
    # In a six-membered ring, the maximum distance is 3; we expect one candidate at distance 3 (position 4)
    candidate_positions = [idx for idx, d in distances.items() if d == 3 and idx not in fused_benzene]
    if len(candidate_positions) == 0:
        return False, "No candidate atom at distance 3 from oxygen (position 4) found."
    elif len(candidate_positions) > 1:
        return False, "Multiple candidate atoms at distance 3 from oxygen found; ambiguous position assignment."
    pos4_atom_idx = candidate_positions[0]
    
    # Step 5: Check that the candidate atom (position 4) carries exactly one substituent that is external to the fused core
    # and that substituent is part of an independent aromatic ring.
    pos4_atom = mol.GetAtomWithIdx(pos4_atom_idx)
    external_substituents = []
    for nbr in pos4_atom.GetNeighbors():
        nbr_idx = nbr.GetIdx()
        if nbr_idx in fused_pyran:
            continue  # Skip atoms that are part of the fused core
        # Ensure the bond is a single bond.
        bond = mol.GetBondBetweenAtoms(pos4_atom_idx, nbr_idx)
        if bond.GetBondTypeAsDouble() != 1:
            continue
        # Check if the neighbor atom is aromatic.
        if not nbr.GetIsAromatic():
            continue
        # Confirm the neighbor is part of an aromatic ring (of at least six atoms) that does not overlap with the fused core.
        for ring in atom_rings:
            if nbr_idx in ring and len(ring) >= 6 and fused_pyran.isdisjoint(ring):
                external_substituents.append(nbr_idx)
                break

    if len(external_substituents) == 0:
        return False, "No external aryl substituent found at the candidate position (position 4)."
    elif len(external_substituents) > 1:
        return False, "Multiple external substituents found at the candidate position; ambiguous assignment."
    
    return True, "Molecule contains a 1-benzopyran core with a unique aryl substituent at position 4; classified as a neoflavonoid."

if __name__ == "__main__":
    # Example test: 4-phenylcoumarin, a prototype neoflavonoid.
    test_smiles = "C12=CC=CC=C1C(=CC(O2)=O)C3=CC=CC=C3"
    result, reason = is_neoflavonoid(test_smiles)
    print("Is neoflavonoid?", result)
    print("Reason:", reason)