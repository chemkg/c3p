"""
Classifies: CHEBI:50128 biflavonoid
"""
"""
Classifies: biflavonoid
Definition:
  A flavonoid oligomer that is obtained by the oxidative coupling of at least two units of aryl‐substituted benzopyran rings
  or its substituted derivatives, resulting in the two ring systems being joined together by a single atom or bond.
  
This approach detects biflavonoids by using two looser SMARTS patterns that look for a benzopyran core:
  - Pattern 1 is a generic aromatic benzopyran (as in many flavones/flavonols)
  - Pattern 2 is a partially saturated benzopyran (as in many flavan/catechin types)
After obtaining substructure matches from both patterns, overlapping matches are merged to ensure that each distinct flavonoid unit is counted only once.
Finally, the connectivity between different flavonoid units is checked to assess proper coupling.
Due to the complexity of flavonoid structures and their many derivatizations, this heuristic may still miss some cases.
"""
from rdkit import Chem

def is_biflavonoid(smiles: str):
    """
    Determines if a molecule is a biflavonoid based on its SMILES string.
    
    The algorithm uses two looser SMARTS patterns to detect flavonoid units:
      - An aromatic benzopyran pattern (generic flavonoid core)
      - A semi-saturated benzopyran pattern (flavan/catechin-like core)
    After finding substructure matches, overlapping matches are merged so that each distinct unit is counted only once.
    Finally, connectivity between at least two distinct units is checked by finding at least one bond that links atoms 
    from different units.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is classified as a biflavonoid, False otherwise.
        str: Explanation of the classification reasoning.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define two looser SMARTS for flavonoid units.
    # Pattern 1: Generic aromatic benzopyran: a benzene ring fused to an oxygen-containing six-membered ring.
    # This should pick up many flavone/flavonol cores.
    aromatic_benzopyran = "c1ccc2Oc(c2c1)"
    # Pattern 2: A semi-saturated variant typical of flavan/catechin units.
    semi_sat_benzopyran = "c1ccc2C(O)Cc2c1"
    
    patt_arom = Chem.MolFromSmarts(aromatic_benzopyran)
    patt_semi = Chem.MolFromSmarts(semi_sat_benzopyran)
    if patt_arom is None and patt_semi is None:
        return False, "Invalid SMARTS patterns"
    
    # Get substructure matches for both patterns.
    matches_arom = mol.GetSubstructMatches(patt_arom)
    matches_semi = mol.GetSubstructMatches(patt_semi)
    
    # Combine the matches (each match is a tuple of atom indices)
    all_matches = list(matches_arom) + list(matches_semi)
    if len(all_matches) < 2:
        return False, "Fewer than two flavonoid-like substructures were found"
    
    # Deduplicate overlapping matches.
    distinct_units = []
    for match in all_matches:
        match_set = set(match)
        found_overlap = False
        for unit in distinct_units:
            overlap = match_set.intersection(unit)
            if len(overlap) >= 0.5 * min(len(match_set), len(unit)):
                found_overlap = True
                unit.update(match_set)  # merge the two matches
                break
        if not found_overlap:
            distinct_units.append(match_set)
    
    if len(distinct_units) < 2:
        return False, "Less than two distinct flavonoid units detected after deduplication"
    
    # Check connectivity: look for at least one bond that connects atoms belonging to two different units.
    connection_found = False
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtomIdx()
        a2 = bond.GetEndAtomIdx()
        units_involved = set()
        for i, unit in enumerate(distinct_units):
            if a1 in unit or a2 in unit:
                units_involved.add(i)
        # If the bond connects atoms from two different distinct units, we have a coupling.
        if len(units_involved) >= 2:
            connection_found = True
            break

    if not connection_found:
        return False, "Two distinct flavonoid units were found, but no coupling bond/atom was detected."
    
    return True, "The molecule contains at least two coupled flavonoid units, consistent with a biflavonoid structure."

# Example usage (for testing):
if __name__ == "__main__":
    # Test with one of the provided SMILES: chamaejasmenin A
    test_smiles = "[H][C@]1(Oc2cc(O)cc(O)c2C(=O)[C@@]1([H])[C@@]1([H])C(=O)c2c(O)cc(O)cc2O[C@]1([H])c1ccc(OC)cc1)c1ccc(OC)cc1"
    result, reason = is_biflavonoid(test_smiles)
    print("Biflavonoid classification:", result)
    print("Reason:", reason)