"""
Classifies: CHEBI:50128 biflavonoid
"""
"""
Classifies: biflavonoid
Definition:
  A flavonoid oligomer that is obtained by the oxidative coupling of at least two units 
  of aryl‐substituted benzopyran rings or its substituted derivatives, resulting in the two 
  ring systems being joined together by a single atom or bond.
  
This approach uses two slightly looser SMARTS patterns (one for an aromatic benzopyran core 
and one for a semi‐saturated variant) to capture a broader range of flavonoid units.
After obtaining substructure matches from both patterns, overlapping matches are merged only 
if they share at least 2 common atoms (rather than a ratio-based threshold). This ensures that 
distinct but coupled units (which may share one atom or bond) are counted separately.
Finally, the connectivity between different flavonoid units is checked.
Due to the enormous diversity in flavonoid derivatizations, this heuristic approach may not be perfect.
"""
from rdkit import Chem

def is_biflavonoid(smiles: str):
    """
    Determines if a molecule is a biflavonoid based on its SMILES string.
    
    The algorithm uses two broadened SMARTS patterns for flavonoid units:
      - An aromatic benzopyran pattern (typical of many flavonol and flavone cores)
      - A semi-saturated benzopyran pattern (typical of many flavan/catechin cores)
      
    Substructure matches are first obtained from both patterns. Then overlapping matches are merged;
    here, matches are merged only when they share at least 2 common atoms, so that coupled flavonoid units 
    (typically sharing only one atom or bond) remain separate.
    Finally, we require that at least two distinct flavonoid units are coupled by at least one bond.
    
    Args:
        smiles (str): SMILES string representing the molecule.
        
    Returns:
        bool: True if the molecule is classified as a biflavonoid, False otherwise.
        str: Explanation of the classification reasoning.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define two broadened SMARTS patterns to capture flavonoid-like units.
    # Pattern 1: Aromatic benzopyran (broadened variant).
    aromatic_benzopyran = "c1ccc2c(c1)oc(c2)"  
    # Pattern 2: Semi-saturated benzopyran (e.g. found in catechins).
    semi_sat_benzopyran = "c1ccc2C(O)Cc2c1"
    
    patt_arom = Chem.MolFromSmarts(aromatic_benzopyran)
    patt_semi = Chem.MolFromSmarts(semi_sat_benzopyran)
    
    if patt_arom is None and patt_semi is None:
        return False, "Invalid SMARTS patterns"
    
    # Get substructure matches from both patterns.
    matches_arom = mol.GetSubstructMatches(patt_arom) if patt_arom is not None else []
    matches_semi = mol.GetSubstructMatches(patt_semi) if patt_semi is not None else []
    
    all_matches = list(matches_arom) + list(matches_semi)
    if len(all_matches) < 2:
        return False, "Fewer than two flavonoid-like substructures were found"
    
    # Deduplicate overlapping matches.
    # Instead of using a ratio, we merge only if the two match sets share at least 2 atoms.
    distinct_units = []
    for match in all_matches:
        match_set = set(match)
        merged = False
        for unit in distinct_units:
            # If the intersection of atom indices is 2 or more, we assume they represent the same core.
            if len(match_set.intersection(unit)) >= 2:
                unit.update(match_set)  # merge sets
                merged = True
                break
        if not merged:
            distinct_units.append(match_set)
            
    if len(distinct_units) < 2:
        return False, "Less than two distinct flavonoid units detected after deduplication"
    
    # Check for connectivity between at least two distinct units:
    # We look for at least one bond that links atoms from two different distinct units.
    connection_found = False
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtomIdx()
        a2 = bond.GetEndAtomIdx()
        units_involved = set()
        for i, unit in enumerate(distinct_units):
            if a1 in unit or a2 in unit:
                units_involved.add(i)
        if len(units_involved) >= 2:
            connection_found = True
            break
    
    if not connection_found:
        return False, "Two distinct flavonoid units were found, but no coupling bond/atom was detected."
    
    return True, "The molecule contains at least two coupled flavonoid units, consistent with a biflavonoid structure."

# Example usage (for testing purposes):
if __name__ == "__main__":
    # Example: chamaejasmenin A SMILES (one of the biflavonoids mentioned)
    test_smiles = "[H][C@]1(Oc2cc(O)cc(O)c2C(=O)[C@@]1([H])[C@@]1([H])C(=O)c2c(O)cc(O)cc2O[C@]1([H])c1ccc(OC)cc1)c1ccc(OC)cc1"
    result, reason = is_biflavonoid(test_smiles)
    print("Biflavonoid classification:", result)
    print("Reason:", reason)