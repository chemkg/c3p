"""
Classifies: CHEBI:35346 11beta-hydroxy steroid
"""
"""
Classifies: 11β-hydroxy steroid
Definition: Any 11-hydroxy steroid in which the hydroxy group at position 11 has beta-configuration.
The classifier attempts to (1) confirm the molecule has a fused-ring core that resembles the steroid nucleus 
(by finding a connected set of rings that is largely carbocyclic and contains ≥17 carbons),
and (2) detect a beta-configured hydroxy group on one of those ring carbons.
Note: This is only an approximate test.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_11beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is an 11β-hydroxy steroid.
    The algorithm:
      1. Parses the SMILES and confirms a minimally steroid-like molecule (≥17 carbons overall).
      2. Identifies the fused ring systems from the molecule’s ring information.
         A steroid nucleus is approximated as a fused (connected) set of rings that is almost entirely carbocyclic
         and contains at least 17 carbon atoms.
      3. Searches for a beta-configured hydroxy group by looking for a chiral carbon with “[C@@H]([O])”
         and then checks that the candidate carbon belongs to the detected steroid core,
         and that the immediate local environment (neighbors in the core) is consistent with a fused-ring.
         
    Args:
      smiles (str): SMILES string of the molecule.
      
    Returns:
      bool: True if the molecule is classified as a 11β-hydroxy steroid, False otherwise.
      str: A reason for the classification.
    """
    # Parse the SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check global carbon count – steroids normally contain many carbons.
    total_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if total_carbons < 17:
        return False, f"Too few carbon atoms ({total_carbons}) to be a steroid"
    
    # Obtain ring information to get all rings (as tuples of atom indices)
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    if not rings:
        return False, "No rings detected; not steroid-like."
    
    # Build fused ring systems – union rings that share at least one atom.
    # We create a list of sets, each representing a connected fused ring system.
    fused_systems = []
    for ring in rings:
        ring_set = set(ring)
        merged = False
        for comp in fused_systems:
            if ring_set & comp:  # if there is any common atom
                comp.update(ring_set)
                merged = True
                break
        if not merged:
            fused_systems.append(ring_set)
    # It is possible that some rings should be merged if they share atoms with more than one component.
    # We repeat the merging process until no further change.
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
    
    # Look for a fused system that qualifies as a steroid core.
    # Our approximate criteria: it contains at least 17 atoms that are carbons and (ideally) very few heteroatoms.
    steroid_core = None
    for comp in fused_systems:
        carbons_in_comp = [idx for idx in comp if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6]
        if len(carbons_in_comp) >= 17:
            # Optionally, check minor contamination by non-carbons.
            non_carbons = len(comp) - len(carbons_in_comp)
            if non_carbons <= 2:  # allow up to 2 heteroatoms in the core
                steroid_core = set(comp)
                break
    if steroid_core is None:
        return False, "No fused ring system with a predominantly carbocyclic core (≥17 carbons) found"
        
    # Now search for a beta-configured hydroxy group.
    # We use a SMARTS pattern for a chiral carbon with '@@' that is attached to an -OH.
    # This pattern is only an approximation for the β configuration.
    beta_oh_smarts = "[C@@H]([O])"
    beta_oh_pat = Chem.MolFromSmarts(beta_oh_smarts)
    if beta_oh_pat is None:
        return False, "Error creating SMARTS pattern"
    
    matches = mol.GetSubstructMatches(beta_oh_pat)
    if not matches:
        return False, "No beta-configured hydroxy group ([C@@H]([O])) detected"
    
    # Now check if any of the matching chiral carbons reside within the steroid core.
    for match in matches:
        carbon_idx = match[0]
        if carbon_idx in steroid_core:
            # Further require that this carbon is bonded to at least two other carbon atoms from the core.
            atom = mol.GetAtomWithIdx(carbon_idx)
            core_neighbors = 0
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in steroid_core and nbr.GetAtomicNum() == 6:
                    core_neighbors += 1
            if core_neighbors >= 2:
                return True, ("Molecule has a fused, predominantly carbocyclic steroid nucleus "
                              "and a beta-configured hydroxy group on a ring carbon (candidate for 11β-hydroxy)")
    
    # If no candidate is found in the steroid core, return False.
    return False, "Found beta-configured hydroxy group(s), but none appear in a suitable steroid nucleus"

# Example usage:
if __name__ == "__main__":
    # Test with one example: cortisol.
    test_smiles = "[H][C@@]12CCC3=CC(=O)CC[C@]3(C)[C@@]1([H])[C@@H](O)C[C@@]1(C)[C@@]2([H])CC[C@]1(O)C(=O)CO"
    result, reason = is_11beta_hydroxy_steroid(test_smiles)
    print(f"Result: {result}\nReason: {reason}")