"""
Classifies: CHEBI:68472 pyrimidine deoxyribonucleoside
"""
"""
Classifies: A deoxyribonucleoside containing a pyrimidine base.
A pyrimidine deoxyribonucleoside contains a five‐membered deoxyribose (furanose) sugar 
that is linked via a glycosidic bond to a pyrimidine base.
"""
from rdkit import Chem

def is_pyrimidine_deoxyribonucleoside(smiles: str):
    """
    Determines if a molecule is a pyrimidine deoxyribonucleoside.
    Heuristic steps:
      1. Detect a pyrimidine base using a SMARTS pattern ("c1cncnc1").
      2. Look for a five-membered ring (furanose) that has exactly one ring oxygen and 4 carbon atoms.
         Then count the exocyclic oxygen substituents on ring carbons (those that are not in the ring and have attached hydrogen(s)).
         In a deoxyribose, we expect exactly two such −OH groups (the 3'-OH and the 5'-CH2OH; the 2' position is deoxy).
      3. Confirm that at least one atom of the sugar ring is directly connected to a nitrogen atom that belongs
         to the pyrimidine base (indicating a glycosidic bond).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a pyrimidine deoxyribonucleoside, False otherwise.
        str: Reason for the classification outcome.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Step 1. Find the pyrimidine base using a SMARTS pattern.
    pyrimidine_smarts = "c1cncnc1"
    pyrimidine_query = Chem.MolFromSmarts(pyrimidine_smarts)
    if not mol.HasSubstructMatch(pyrimidine_query):
        return False, "Pyrimidine base not found"
    pyrimidine_matches = mol.GetSubstructMatches(pyrimidine_query)
    # Create a set of all atom indices in any pyrimidine match.
    pyrimidine_atoms = set(idx for match in pyrimidine_matches for idx in match)

    # Step 2. Identify candidate deoxyribose sugar rings.
    ring_info = mol.GetRingInfo()
    candidate_sugar = None
    for ring in ring_info.AtomRings():
        if len(ring) != 5:
            continue  # Skip non-5-membered rings.
        # Count heavy atoms in the ring by their atomic symbol.
        n_ring_oxy = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetSymbol() == "O")
        n_ring_carbons = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetSymbol() == "C")
        # For a furanose ring, we expect exactly one oxygen and four carbons.
        if n_ring_oxy != 1 or n_ring_carbons != 4:
            continue

        # Count exocyclic oxygen substituents on ring carbons.
        exocylic_OH = 0
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetSymbol() != "C":
                continue  # Only check for oxygens on carbons.
            for nbr in atom.GetNeighbors():
                # Consider only neighbors that are not part of the ring.
                if nbr.GetIdx() in ring:
                    continue
                if nbr.GetSymbol() == "O":
                    # Check if this oxygen atom has at least one hydrogen attached.
                    # Using GetTotalNumHs() returns the number of (implicit+explicit) H.
                    if nbr.GetTotalNumHs() > 0:
                        exocylic_OH += 1
        # We expect exactly 2 exocyclic OH groups (the deoxyribose lacks the 2'-OH).
        if exocylic_OH != 2:
            continue

        # Step 3. Verify glycosidic linkage:
        # Check that one of the ring atoms is directly bonded to a nitrogen that belongs to the pyrimidine.
        sugar_linked = False
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in pyrimidine_atoms and nbr.GetSymbol() == "N":
                    sugar_linked = True
                    break
            if sugar_linked:
                break

        if not sugar_linked:
            continue  # This candidate ring is not linked to the base.
            
        # If we reach here, we have found a candidate deoxyribose sugar moiety.
        candidate_sugar = ring
        break

    if candidate_sugar is None:
        return False, "Deoxyribose sugar moiety not found"

    return True, "Molecule contains a pyrimidine base linked to a deoxyribose sugar via a glycosidic bond"

# For testing purposes, you can call the function with one or several provided SMILES.
if __name__ == "__main__":
    # Example: 2'-deoxyuridine
    test_smiles = "OC[C@H]1O[C@H](C[C@@H]1O)n1ccc(=O)[nH]c1=O"
    result, reason = is_pyrimidine_deoxyribonucleoside(test_smiles)
    print(result, reason)