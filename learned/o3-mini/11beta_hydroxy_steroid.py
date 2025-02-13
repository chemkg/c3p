"""
Classifies: CHEBI:35346 11beta-hydroxy steroid
"""
"""
Classifies: 11β-hydroxy steroid
Definition: A fused polycyclic (predominantly carbocyclic) steroid nucleus (with at least 3 rings and ≥17 carbons)
with a beta-configured hydroxy group. The beta-configured hydroxy is approximated by detecting a chiral center
that is described by the SMARTS pattern "[C@@H]([O])", but only if that chiral carbon is part of a steroid-like core.
Note: This is still only an approximate test.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_11beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is an 11β-hydroxy steroid.
    The algorithm:
      1. Parses the SMILES and checks that there are enough carbon atoms.
      2. Extracts all rings and builds fused (connected) ring systems.
      3. Among fused systems, selects one that is predominantly carbocyclic,
         contains at least 17 carbon atoms AND is made up of at least 3 rings.
      4. Searches for beta-configured hydroxy groups using the SMARTS pattern "[C@@H]([O])".
      5. For any candidate chiral carbon, confirms that
           - it belongs to the selected steroid core,
           - it is in a ring (and one of its rings is typical in size for steroids, i.e. 5 or 6 members),
           - and it has at least two carbon neighbors that are also in the core.
      
    Args:
      smiles (str): SMILES string of the molecule.
      
    Returns:
      bool: True if the molecule qualifies as an 11β-hydroxy steroid, False otherwise.
      str: A reason text for the classification.
    """
    # Parse the SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Global check: count total carbons; steroids normally have many carbons.
    total_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if total_carbons < 17:
        return False, f"Too few carbon atoms ({total_carbons}) to be a steroid"

    # Obtain ring information: list of tuples of atom indices for every ring.
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    if not rings:
        return False, "No rings detected; not steroid-like."
    
    # Build fused ring systems: each fused system is a set containing all atom indices in connected rings.
    fused_systems = []
    for ring in rings:
        ring_set = set(ring)
        merged = False
        for comp in fused_systems:
            if ring_set & comp:  # if overlap exists
                comp.update(ring_set)
                merged = True
                break
        if not merged:
            fused_systems.append(ring_set)
    # It is possible that some systems need further merging.
    changed = True
    while changed:
        changed = False
        new_systems = []
        while fused_systems:
            comp = fused_systems.pop(0)
            merged_found = False
            for i, other in enumerate(fused_systems):
                if comp & other:
                    fused_systems[i] = comp | other
                    merged_found = True
                    changed = True
                    break
            if not merged_found:
                new_systems.append(comp)
        fused_systems = new_systems[:]
    
    # From the fused systems, look for one that could represent a steroid nucleus.
    steroid_core = None
    for comp in fused_systems:
        # Count how many carbons are present in the fused system.
        carbons_in_comp = [idx for idx in comp if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6]
        if len(carbons_in_comp) < 17:
            continue
        # Additionally, require the system to incorporate at least 3 rings.
        ring_count = 0
        for ring in rings:
            if set(ring).issubset(comp):
                ring_count += 1
        if ring_count < 3:
            continue
        # Optionally allow very few heteroatoms in the core.
        non_carbons = len(comp) - len(carbons_in_comp)
        if non_carbons > 2:
            continue
        steroid_core = set(comp)
        break

    if steroid_core is None:
        return False, "No fused ring system with a predominantly carbocyclic core (≥17 carbons and ≥3 rings) found"

    # Prepare the SMARTS for beta-configured hydroxy: chiral carbon with an -OH.
    beta_oh_smarts = "[C@@H]([O])"
    beta_oh_pat = Chem.MolFromSmarts(beta_oh_smarts)
    if beta_oh_pat is None:
        return False, "Error creating SMARTS pattern for beta-configured hydroxyl"

    matches = mol.GetSubstructMatches(beta_oh_pat)
    if not matches:
        return False, "No beta-configured hydroxy group ([C@@H]([O])) detected"

    # Check if any matching chiral carbon appears in the steroid core.
    # Additionally, require that the candidate is in a ring (of size 5 or 6) and has ≥2 core carbon neighbors.
    # This should reduce false positives where the beta-OH occurs in non-steroidal contexts.
    for match in matches:
        carbon_idx = match[0]
        if carbon_idx not in steroid_core:
            continue
        atom = mol.GetAtomWithIdx(carbon_idx)
        if not atom.IsInRing():
            continue

        # Check if the candidate belongs to a ring of plausible steroid size (5 or 6 members)
        candidate_in_valid_ring = False
        for ring in rings:
            if carbon_idx in ring and (len(ring) == 5 or len(ring) == 6):
                # further require that the entire ring is within the steroid core
                if set(ring).issubset(steroid_core):
                    candidate_in_valid_ring = True
                    break
        if not candidate_in_valid_ring:
            continue

        # Count the number of neighbors from the steroid core that are carbons.
        core_neighbors = 0
        for nbr in atom.GetNeighbors():
            if nbr.GetIdx() in steroid_core and nbr.GetAtomicNum() == 6:
                core_neighbors += 1
        if core_neighbors >= 2:
            return True, ("Molecule has a fused, predominantly carbocyclic steroid nucleus (≥3 rings, ≥17 carbons) "
                          "and a beta-configured hydroxy group on a ring carbon (candidate for 11β-hydroxy).")
    
    return False, "Found beta-configured hydroxy group(s), but none appear in a suitable steroid nucleus"

# Example usage:
if __name__ == "__main__":
    # Test with one example (cortisol):
    test_smiles = "[H][C@@]12CCC3=CC(=O)CC[C@]3(C)[C@@]1([H])[C@@H](O)C[C@@]1(C)[C@@]2([H])CC[C@]1(O)C(=O)CO"
    result, reason = is_11beta_hydroxy_steroid(test_smiles)
    print(f"Result: {result}\nReason: {reason}")