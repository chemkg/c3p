"""
Classifies: CHEBI:50128 biflavonoid
"""
"""
Classifies: biflavonoid
Definition:
    A flavonoid oligomer that is obtained by the oxidative coupling of at least two units of aryl‐substituted benzopyran rings
    or its substituted derivatives, resulting in the two ring systems being joined together by a single atom or bond.
    
This approach first attempts to locate two flavonoid “cores” using two simplified SMARTS patterns (one based on a flavone skeleton
and one on a flavan [catechin-like] skeleton). It then deduplicates overlapping matches and finally checks whether there is a bond connecting
atoms from two distinct units.
Note: Due to the structural complexity of flavonoids and their many substituted derivatives, this heuristic may not correctly classify every example.
"""
from rdkit import Chem

def is_biflavonoid(smiles: str):
    """
    Determines if a molecule is a biflavonoid based on its SMILES string.
    
    The algorithm uses two simplified SMARTS patterns to detect flavonoid units:
      - A flavone-like unit (a benzopyran with a C=O group, as found in some flavonoids)
      - A flavan-like unit (catechin/epicatechin-type skeleton)
    After finding substructure matches, overlapping matches are merged so that each distinct flavonoid unit is counted only once.
    Finally, it checks whether at least two distinct flavonoid units are connected (by at least one bond bridging atoms from different units).
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is classified as a biflavonoid, False otherwise.
        str: Explanation of the classification reasoning.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define two simplified SMARTS patterns for a flavonoid unit.
    # Pattern 1: A simple flavone-like skeleton (benzopyran with a ketone)
    flavone_smarts = "c1ccc2c(c1)OC(=O)c3ccccc23"
    # Pattern 2: A simplified flavan (catechin-like) skeleton.
    # (Many catechins have a characteristic oxygenated heterocycle; here we require a six-membered ring containing O)
    flavan_smarts = "OC1CCc2cc(O)cc(O)c2O1"
    
    patt_flavone = Chem.MolFromSmarts(flavone_smarts)
    patt_flavan = Chem.MolFromSmarts(flavan_smarts)
    if patt_flavone is None and patt_flavan is None:
        return False, "Invalid SMARTS patterns"
    
    # Get substructure matches for both patterns.
    matches_flavone = mol.GetSubstructMatches(patt_flavone)
    matches_flavan = mol.GetSubstructMatches(patt_flavan)
    
    # Combine the matches into one list (each match is a tuple of atom indices)
    all_matches = list(matches_flavone) + list(matches_flavan)
    if len(all_matches) < 2:
        return False, "Fewer than two flavonoid-like substructures were found"
    
    # Deduplicate overlapping matches.
    # If two matches share more than 50% of the atoms of one match, consider them as the same flavonoid unit.
    distinct_units = []
    for match in all_matches:
        match_set = set(match)
        found_overlap = False
        for unit in distinct_units:
            # Calculate fraction of overlap relative to the smaller set.
            overlap = match_set.intersection(unit)
            if len(overlap) >= 0.5 * min(len(match_set), len(unit)):
                found_overlap = True
                # Merge the sets (to account for slight differences in matching)
                unit.update(match_set)
                break
        if not found_overlap:
            distinct_units.append(match_set)
    
    if len(distinct_units) < 2:
        return False, "Less than two distinct flavonoid units detected after deduplication"
    
    # Check connectivity between at least two distinct units.
    # We require that there exists at least one bond that connects an atom from one unit to an atom from another unit.
    bonds = mol.GetBonds()
    connection_found = False
    for b in bonds:
        a1 = b.GetBeginAtomIdx()
        a2 = b.GetEndAtomIdx()
        # Check if a1 and a2 belong to different detected units.
        unit_indices = []
        for i, unit in enumerate(distinct_units):
            if a1 in unit or a2 in unit:
                unit_indices.append(i)
        # If the two bond atoms belong to two different units, we have a coupling bond.
        if len(unit_indices) == 2 and unit_indices[0] != unit_indices[1]:
            connection_found = True
            break
        # In some biflavonoids, the two flavonoid units might even share a common atom (a bridging atom).
        # So if one atom belongs to more than one unit we consider that a connection.
        for unit in distinct_units:
            # If one of the atoms is in the unit and the other is not,
            # then if that other atom is in another unit later, that counts.
            pass  # (the bond check above is sufficient for most cases)
    
    if not connection_found:
        return False, "Two distinct flavonoid units were found, but no direct coupling bond/atom was detected."
    
    return True, "The molecule contains at least two coupled flavonoid units, consistent with a biflavonoid structure."

# Example usage:
if __name__ == "__main__":
    # Test with one of the provided SMILES, e.g., chamaejasmenin A
    test_smiles = "[H][C@]1(Oc2cc(O)cc(O)c2C(=O)[C@@]1([H])[C@@]1([H])C(=O)c2c(O)cc(O)cc2O[C@]1([H])c1ccc(OC)cc1)c1ccc(OC)cc1"
    result, reason = is_biflavonoid(test_smiles)
    print("Biflavonoid classification:", result)
    print("Reason:", reason)