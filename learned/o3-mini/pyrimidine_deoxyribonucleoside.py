"""
Classifies: CHEBI:68472 pyrimidine deoxyribonucleoside
"""
"""
Classifies: A deoxyribonucleoside containing a pyrimidine base.
A pyrimidine deoxyribonucleoside is defined as a molecule that has a pyrimidine base
(canonical SMARTS used here) linked via a glycosidic bond to a five‐membered deoxyribose sugar.
We also reject molecules containing phosphorus (i.e. nucleotides) and enforce that the sugar
ring is “deoxy” by checking that its total oxygen count (ring oxygen plus oxygen substituents)
is either 3 or 4 (but not 5 or more, as that is typical of ribose).
"""

from rdkit import Chem

def is_pyrimidine_deoxyribonucleoside(smiles: str):
    """
    Determines if a molecule is a pyrimidine deoxyribonucleoside.
    
    Heuristic steps:
      1. Reject molecules containing phosphorus (likely nucleotides).
      2. Identify a pyrimidine base by a simple SMARTS pattern.
      3. Identify candidate five-membered rings (sugar ring) that have exactly one ring oxygen and four ring carbons.
         Then for each candidate, count the exocyclic oxygen atoms attached to ring carbons.
         In a canonical deoxyribose the total oxygen count (ring oxygen + exocyclic oxygens) is 4;
         however, we allow a total count of 3 if one -OH (typically the 3'-OH) is replaced by another group.
         We reject candidates with a total oxygen count of 5 or more.
      4. Confirm that at least one of the ring carbons is directly bonded to a nitrogen atom that is part
         of the pyrimidine base (indicative of a glycosidic linkage).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a pyrimidine deoxyribonucleoside, False otherwise.
        str: A reason explaining the outcome.
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Step 1: Exclude molecules with phosphorus (P) – likely nucleotides not nucleosides.
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == "P":
            return False, "Phosphate group detected – molecule is likely a nucleotide"
    
    # Step 2: Identify the pyrimidine base.
    # Use a simple pyrimidine SMARTS pattern.
    pyrimidine_smarts = "c1cncnc1"
    pyrimidine_query = Chem.MolFromSmarts(pyrimidine_smarts)
    if not mol.HasSubstructMatch(pyrimidine_query):
        return False, "Pyrimidine base not found"
    pyrimidine_matches = mol.GetSubstructMatches(pyrimidine_query)
    # Aggregate all atom indices from the pyrimidine part
    pyrimidine_atoms = set(idx for match in pyrimidine_matches for idx in match)
    
    # Step 3: Search for candidate sugar rings.
    ring_info = mol.GetRingInfo()
    candidate_sugar = None
    sugar_total_o_count = None
    for ring in ring_info.AtomRings():
        if len(ring) != 5:
            continue  # must be five-membered ring
            
        # Count ring oxygens and carbons.
        ring_oxygens = [idx for idx in ring if mol.GetAtomWithIdx(idx).GetSymbol() == "O"]
        ring_carbons = [idx for idx in ring if mol.GetAtomWithIdx(idx).GetSymbol() == "C"]
        if len(ring_oxygens) != 1 or len(ring_carbons) != 4:
            continue
            
        # Count exocyclic oxygen atoms attached to ring carbons.
        exo_oxygen_count = 0
        for idx in ring_carbons:
            atom = mol.GetAtomWithIdx(idx)
            for nbr in atom.GetNeighbors():
                # Only count oxygen atoms that are not part of the ring.
                if nbr.GetIdx() in ring:
                    continue
                if nbr.GetSymbol() == "O":
                    # We require that the oxygen is not double bonded (likely an -OH or -O- in an ether)
                    # and that it has at least one hydrogen (an -OH group typically).
                    if nbr.GetTotalNumHs() > 0:
                        exo_oxygen_count += 1
        total_o = len(ring_oxygens) + exo_oxygen_count  # total oxygen count in candidate sugar
        
        # In ribose (non-deoxy) total oxygens would be 5 or more.
        # In deoxyribose we expect a total of 4 in the canonical case; however, modifications (e.g. thio substitution)
        # might reduce the count to 3. We accept total oxygen count of 3 or 4.
        if total_o not in (3, 4):
            continue
            
        # Step 4: Confirm that the sugar is linked to the pyrimidine via a glycosidic bond.
        # Look for any ring atom (preferably a carbon) that is directly bonded to a nitrogen belonging to the pyrimidine.
        sugar_linked = False
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            # We usually check only carbons, as the glycosidic link is almost always from a carbon.
            if atom.GetSymbol() != "C":
                continue
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in pyrimidine_atoms and nbr.GetSymbol() == "N":
                    sugar_linked = True
                    break
            if sugar_linked:
                break
                
        if not sugar_linked:
            continue  # this ring is not glycosidically linked
        
        # This candidate sugar ring qualifies.
        candidate_sugar = ring
        sugar_total_o_count = total_o
        break

    if candidate_sugar is None:
        return False, "Deoxyribose sugar moiety not found (or sugar oxygen count indicates ribose)"
        
    return True, f"Molecule contains a pyrimidine base and a deoxyribose sugar (total sugar O count={sugar_total_o_count}) linked by a glycosidic bond"

# For testing the function you can run this module directly.
if __name__ == "__main__":
    # Example: 2'-deoxyuridine
    test_smiles = "OC[C@H]1O[C@H](C[C@@H]1O)n1ccc(=O)[nH]c1=O"
    result, reason = is_pyrimidine_deoxyribonucleoside(test_smiles)
    print(result, reason)