"""
Classifies: CHEBI:35346 11beta-hydroxy steroid
"""
"""
Classifies: 11β-hydroxy steroid
Definition: A fused, polycyclic, predominantly carbocyclic steroid nucleus – ideally the classic tetracyclic (cyclopentanoperhydrophenanthrene) system –
which consists of four fused rings (typically three six‐membered and one five‐membered, or an aromatic variant) 
with at least 17 carbons in the fused system, along with at least one beta‑configured hydroxy group (detected by a chiral center with –OH in beta configuration).
Note: This is still an approximate test.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_11beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is an 11β-hydroxy steroid.
    
    The algorithm is as follows:
      1. Parse the SMILES and ensure the molecule has enough carbons (≥17).
      2. Extract ring data and build fused ring systems. A fused system is a set of rings sharing at least one atom.
      3. Among the fused systems, look for one that is most likely the steroid nucleus – ideally a tetracyclic core.
         We require that one fused system contains exactly four rings whose union has at least 17 carbon atoms.
         Moreover, the rings should represent either the classic 3 six‐membered + 1 five‐membered steroid skeleton,
         or – in cases such as aromatic estrogens – 4 six‐membered rings.
      4. Use a SMARTS pattern "[C@@H][OX2H]" to detect a beta‑configured hydroxy group and then check that
         the candidate chiral carbon is part of the steroid core (in a ring of size 5 or 6 that is completely within said core)
         and has at least two neighbors (that are carbons) in the core.
    
    Args:
      smiles (str): SMILES string of the molecule.
      
    Returns:
      bool: True if the molecule qualifies as an 11β-hydroxy steroid, False otherwise.
      str: A reason text for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Global carbon count check.
    total_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if total_carbons < 17:
        return False, f"Too few carbons ({total_carbons}) to be a steroid"
    
    # Get ring information.
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    if not rings:
        return False, "No rings found; not steroid-like."
    
    # Build fused ring systems (each is a set of atom indices from connected rings).
    fused_systems = []
    for ring in rings:
        ring_set = set(ring)
        merged = False
        for comp in fused_systems:
            if ring_set & comp:  # if share any atoms, merge
                comp.update(ring_set)
                merged = True
                break
        if not merged:
            fused_systems.append(ring_set)
    # Further merge any overlapping systems.
    changed = True
    while changed:
        changed = False
        new_systems = []
        while fused_systems:
            comp = fused_systems.pop(0)
            merge_found = False
            for i, other in enumerate(fused_systems):
                if comp & other:
                    fused_systems[i] = comp | other
                    merge_found = True
                    changed = True
                    break
            if not merge_found:
                new_systems.append(comp)
        fused_systems = new_systems[:]
    
    # Look among fused systems for a steroid nucleus candidate.
    steroid_core = None
    for comp in fused_systems:
        # Identify rings fully inside this fused system.
        core_rings = [ring for ring in rings if set(ring).issubset(comp)]
        if len(core_rings) != 4:
            continue  # prefer tetracyclic systems
        # Count carbons in the fused system.
        comp_carbons = [idx for idx in comp if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6]
        if len(comp_carbons) < 17:
            continue
        # Determine ring sizes within the core.
        ring_sizes = [len(ring) for ring in core_rings]
        count5 = ring_sizes.count(5)
        count6 = ring_sizes.count(6)
        # Accept if either the classical steroid (1 five-membered and 3 six-membered) or a variant (4 six-membered).
        if (count5 == 1 and count6 == 3) or (count5 == 0 and count6 == 4):
            steroid_core = comp
            break
    
    # Fallback: if no tetracyclic core is found, try any fused system with at least 3 rings.
    if steroid_core is None:
        for comp in fused_systems:
            core_rings = [ring for ring in rings if set(ring).issubset(comp)]
            if len(core_rings) < 3:
                continue
            comp_carbons = [idx for idx in comp if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6]
            if len(comp_carbons) < 17:
                continue
            steroid_core = comp
            break

    if steroid_core is None:
        return False, "No suitable fused ring system (steroid nucleus) found with enough carbons and rings."
    
    # Prepare the SMARTS for beta-configured hydroxyl (the beta-OH is approximated with a chiral carbon [C@@H] directly bonded to –OH).
    beta_oh_smarts = "[C@@H][OX2H]"
    beta_oh_pat = Chem.MolFromSmarts(beta_oh_smarts)
    if beta_oh_pat is None:
        return False, "Error compiling SMARTS for beta-configured hydroxyl."
    
    matches = mol.GetSubstructMatches(beta_oh_pat)
    if not matches:
        return False, "No beta‑configured –OH (pattern [C@@H][OX2H]) detected."
    
    # For each candidate beta‑OH, verify that the chiral carbon is in a valid position in the steroid core.
    for match in matches:
        candidate_idx = match[0]  # the chiral carbon index
        if candidate_idx not in steroid_core:
            continue  # candidate not in the steroid core
        candidate_atom = mol.GetAtomWithIdx(candidate_idx)
        if not candidate_atom.IsInRing():
            continue  # must be part of a ring
        # Check that the candidate is part of a ring which is fully inside the steroid core and is of plausible size.
        found_valid_ring = False
        for ring in rings:
            if candidate_idx in ring and (len(ring) == 5 or len(ring) == 6):
                if set(ring).issubset(steroid_core):
                    found_valid_ring = True
                    break
        if not found_valid_ring:
            continue
        
        # Check that the candidate carbon has at least two neighboring carbons from the core.
        core_neighbor_count = 0
        for nbr in candidate_atom.GetNeighbors():
            if nbr.GetIdx() in steroid_core and nbr.GetAtomicNum() == 6:
                core_neighbor_count += 1
        if core_neighbor_count >= 2:
            return True, ("Molecule has a tetracyclic (or near‐tetracyclic) steroid nucleus "
                          "and a beta‑configured –OH group "
                          "located in a valid position within that nucleus.")
    
    return False, "Found beta‑configured –OH group(s) but none in a suitable steroid nucleus context."

# Example usage:
if __name__ == "__main__":
    # Here we test the function with one example (cortisol) from the list.
    test_smiles = "[H][C@@]12CCC3=CC(=O)CC[C@]3(C)[C@@]1([H])[C@@H](O)C[C@@]1(C)[C@@]2([H])CC[C@]1(O)C(=O)CO"
    result, reason = is_11beta_hydroxy_steroid(test_smiles)
    print(f"Result: {result}\nReason: {reason}")