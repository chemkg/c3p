"""
Classifies: CHEBI:50128 biflavonoid
"""
"""
Classifies: biflavonoid
Definition:
  A flavonoid oligomer that is obtained by the oxidative coupling of at least two units 
  of aryl‐substituted benzopyran rings or its substituted derivatives, resulting in the two 
  ring systems being joined together by a single atom or bond.

This approach uses two broadened SMARTS patterns for flavonoid units:
  - A flavonol-type (chromen-4-one) core
  - A flavan-type (dihydroflavonoid) core
  
Substructure matches from these patterns are collected and then “deduplicated” only in cases 
when one match is essentially contained within another. This preserves two distinct units 
even when they share a coupling bond or one atom. Then, we require that at least one bond 
connects atoms from two different (deduplicated) units.
"""

from rdkit import Chem

def is_biflavonoid(smiles: str):
    """
    Determines if a molecule is a biflavonoid based on its SMILES string.
    
    The algorithm uses multiple SMARTS patterns to detect flavonoid-like cores:
      - A flavonol core (chromen-4-one skeleton): example pattern "O=c1c(O)cc(O)c2occc12"
      - A flavan (dihydroflavonoid) core: example pattern "c1ccc2C(O)Cc2c1"
      
    Substructure matches from all SMARTS patterns are gathered. Overlapping matches are merged 
    only if one substructure is completely contained in (or contains) the other. This avoids erroneously 
    merging two coupled units (which might share only one atom or bond). A final connectivity check 
    requires that at least one bond connects atoms from two distinct flavonoid-like substructures.
    
    Args:
        smiles (str): SMILES string representing the molecule.
        
    Returns:
        bool: True if the molecule is classified as a biflavonoid, False otherwise.
        str: Explanation of the classification reasoning.
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a list of broadened SMARTS patterns for flavonoid-like cores.
    # Pattern 1: Flavonol/chromen-4-one (common in many flavonoids)
    smarts_core1 = "O=c1c(O)cc(O)c2occc12"
    # Pattern 2: Flavan-3-ol (dihydroflavonoid core found in catechin units)
    smarts_core2 = "c1ccc2C(O)Cc2c1"
    # (Optionally, one can add additional patterns to capture more structural variability.)
    
    patterns = []
    for s in [smarts_core1, smarts_core2]:
        patt = Chem.MolFromSmarts(s)
        if patt is not None:
            patterns.append(patt)
            
    if not patterns:
        return False, "No valid SMARTS patterns defined"
    
    # Obtain substructure matches from all patterns.
    all_matches = []
    for patt in patterns:
        ms = mol.GetSubstructMatches(patt)
        if ms:
            all_matches.extend(ms)
    
    if len(all_matches) < 2:
        return False, "Fewer than two flavonoid-like substructures were found"
    
    # Deduplicate overlapping matches. Instead of merging when two matches share at least 2 atoms,
    # we merge only when one match is entirely (or almost entirely) contained in another.
    distinct_units = []
    for match in all_matches:
        match_set = set(match)
        merged = False
        for i, unit in enumerate(distinct_units):
            # If one set is a subset of the other, assume they represent the same core.
            if match_set.issubset(unit) or unit.issubset(match_set):
                distinct_units[i] = unit.union(match_set)
                merged = True
                break
        if not merged:
            distinct_units.append(match_set)
    
    if len(distinct_units) < 2:
        return False, "Less than two distinct flavonoid units detected after deduplication"
    
    # Check for connectivity between at least two distinct units:
    # This means that at least one bond in the molecule should connect atoms belonging to different units.
    connection_found = False
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtomIdx()
        a2 = bond.GetEndAtomIdx()
        units_found = set()
        for idx, unit in enumerate(distinct_units):
            if a1 in unit or a2 in unit:
                units_found.add(idx)
        if len(units_found) >= 2:
            connection_found = True
            break
            
    if not connection_found:
        return False, "Two distinct flavonoid units were found, but no coupling bond/atom was detected"
    
    return True, "The molecule contains at least two coupled flavonoid units, consistent with a biflavonoid structure"

# Example testing code (if run as a script, uncomment the lines below):
# if __name__ == "__main__":
#     test_smiles = "[H][C@]1(Oc2cc(O)cc(O)c2C(=O)[C@@]1([H])[C@@]1([H])C(=O)c2c(O)cc(O)cc2O[C@]1([H])c1ccc(OC)cc1)c1ccc(OC)cc1"
#     result, reason = is_biflavonoid(test_smiles)
#     print("Biflavonoid classification:", result)
#     print("Reason:", reason)