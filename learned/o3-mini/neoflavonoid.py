"""
Classifies: CHEBI:71971 neoflavonoid
"""
Below is one possible improved solution. In our reasoning we concluded that simply “counting” any external aromatic neighbor of any non‐fused atom in the pyran ring was too unspecific. Instead we try to mimic the “numbering” of a 1–benzopyran. In a proper 1–benzopyran the only heteroatom (oxygen) is at position 1 so that the atom at position 4 (three bonds away in the six‐membered ring) is the only “free” ring vertex that may carry an external aryl substituent. Our new approach is therefore to

  (1) find all six–membered rings with exactly one oxygen (candidate “pyran” rings),
    (2) find benzene rings (fully aromatic six–membered rings),
    (3) choose one candidate pyran ring which is fused (shares at least two atoms) with a benzene ring so that we have a candidate benzopyran core,
    (4) identify the unique oxygen in the candidate pyran, build a connectivity graph on the atoms in the pyran ring, and compute the shortest–path distances.
    (5) select (if unique) the atom at distance 3 from the oxygen and not belonging to the fused benzene portion (this candidate is “position 4”).
    (6) then check that exactly one substituent (via a single bond) is attached to that candidate and that the substituent is part of an aromatic ring that is entirely separate from the candidate core.
If any step fails or if more than one candidate is found we return False with the appropriate message.

Below is the complete Python code in one block:

------------------------------------------------

#!/usr/bin/env python
"""
Classifies: neoflavonoid
Definition: A neoflavonoid is any 1‐benzopyran with an aryl substituent at position 4.
Revised heuristic:
  1. Identify candidate pyran rings: 6‐membered rings containing exactly one oxygen.
  2. Identify candidate benzene rings: aromatic 6‐membered rings.
  3. Find a pyran ring fused (sharing at least two atoms) with a benzene ring.
  4. In that candidate pyran ring, locate the unique oxygen and compute the distances (along bonds within the ring).
     The candidate “position 4” is the atom at distance 3 (the maximum distance in a 6–membered ring) from the oxygen,
     provided it is not also part of the fused benzene ring.
  5. Check that that candidate atom has exactly one substituent (via a single bond) outside of the benzopyran core,
     and that the substituent is part of an independent aromatic ring.
If any step fails or is ambiguous the classification returns False.
"""

from rdkit import Chem

def is_neoflavonoid(smiles: str):
    """
    Determines if a molecule is a neoflavonoid based on its SMILES string.
    A neoflavonoid is defined as any 1-benzopyran with an aryl substituent at position 4.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a neoflavonoid, False otherwise.
        str: Explanation for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    try:
        Chem.SanitizeMol(mol)
    except Exception as e:
        return False, f"Sanitization error: {str(e)}"
    
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    
    # Step 1: Find candidate pyran rings: rings with 6 atoms and exactly one oxygen.
    pyran_rings = []
    for ring in atom_rings:
        if len(ring) == 6:
            oxy_count = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8)
            if oxy_count == 1:
                pyran_rings.append(set(ring))
                
    if not pyran_rings:
        return False, "No six-membered oxygen-containing (pyran) ring found."
    
    # Step 2: Find candidate benzene rings: aromatic rings with 6 atoms.
    benzene_rings = []
    for ring in atom_rings:
        if len(ring) == 6:
            # require all atoms aromatic
            if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
                benzene_rings.append(set(ring))
    
    if not benzene_rings:
        return False, "No aromatic benzene ring found for potential fusion."
    
    # Step 3: Identify a pyran ring that is fused with a benzene ring (sharing at least 2 atoms)
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
        return False, "No fused benzopyran (pyran fused to benzene) core found."
    
    # Step 4: In the candidate pyran ring, identify the unique oxygen atom.
    oxygen_atoms = [idx for idx in fused_pyran if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8]
    if len(oxygen_atoms) != 1:
        return False, "Candidate pyran ring does not contain exactly one oxygen."
    oxy_idx = oxygen_atoms[0]
    
    # Build a simple connectivity graph (only bonds among atoms in fused_pyran)
    connectivity = { idx: set() for idx in fused_pyran }
    for idx in fused_pyran:
        atom = mol.GetAtomWithIdx(idx)
        for nbr in atom.GetNeighbors():
            nbr_idx = nbr.GetIdx()
            if nbr_idx in fused_pyran:
                connectivity[idx].add(nbr_idx)
    
    # Perform a breadth-first search from the oxygen to compute distances
    distances = { oxy_idx: 0 }
    queue = [oxy_idx]
    while queue:
        current = queue.pop(0)
        for nbr in connectivity[current]:
            if nbr not in distances:
                distances[nbr] = distances[current] + 1
                queue.append(nbr)
    # In a six-membered ring the maximum bond distance from any vertex is 3.
    candidate_positions = [idx for idx, dist in distances.items() 
                           if dist == 3 and idx not in fused_benzene]
    if len(candidate_positions) == 0:
        return False, "No candidate position at distance 3 from oxygen (likely position 4) found."
    elif len(candidate_positions) > 1:
        return False, "Multiple candidate positions at distance 3 from oxygen found; ambiguous assignment."
    
    pos4_atom_idx = candidate_positions[0]
    
    # Step 5: Look at neighbors of the candidate atom that lie outside the entire fused (benzopyran) core.
    external_substituents = []
    pos4_atom = mol.GetAtomWithIdx(pos4_atom_idx)
    for nbr in pos4_atom.GetNeighbors():
        nbr_idx = nbr.GetIdx()
        if nbr_idx in fused_pyran:
            continue  # skip atoms in the core
        # We want a substituent that is directly attached (single bond) and is part of an aromatic ring.
        if not nbr.GetIsAromatic():
            continue
        # Check if the neighbor is in an aromatic ring (of size>=6) that is entirely outside the fused core.
        in_external_ring = False
        for ring in atom_rings:
            if nbr_idx in ring and len(ring) >= 6:
                if fused_pyran.intersection(ring) == set():
                    in_external_ring = True
                    break
        if in_external_ring:
            external_substituents.append(nbr_idx)
            
    if len(external_substituents) == 0:
        return False, "No aryl substituent detected at candidate position (likely position 4)."
    elif len(external_substituents) > 1:
        return False, "Multiple substituents at candidate position detected; ambiguous neoflavonoid assignment."
    
    return True, ("Molecule contains a 1-benzopyran core with a unique aryl substituent attached at the candidate "
                   "position (position 4); pattern consistent with a neoflavonoid.")

# Example usage:
if __name__ == "__main__":
    # test with a known benchmark: 4-phenylcoumarin (neoflavonoid prototype)
    test_smiles = "C12=CC=CC=C1C(=CC(O2)=O)C3=CC=CC=C3"
    result, reason = is_neoflavonoid(test_smiles)
    print("Is neoflavonoid?", result)
    print("Reason:", reason)

------------------------------------------------

Notes and caveats:
• This code uses a heuristic that first finds a candidate fused bipartite system (pyran fused to benzene) and then "numbers" the pyran ring by finding the atom three bonds away from the only oxygen.
• Only if that candidate atom bears a single substituent leading into a fully independent aromatic ring (of at least six atoms) do we classify the molecule as a neoflavonoid.
• Molecules that fall outside this strict interpretation, or that contain additional fused rings or ambiguous substituents, will not be classified as neoflavonoids.

This approach (while still heuristic) should reduce both the false positives and negatives reported previously.