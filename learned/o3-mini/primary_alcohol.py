"""
Classifies: CHEBI:15734 primary alcohol
"""
"""
Classifies: Primary Alcohol

Definition:
  A primary alcohol is a compound in which at least one hydroxyl (-OH) group is attached
  to a saturated (sp3) carbon that is “primary” in the sense that it is either 
  (a) a methanol unit (–CH3–OH) or 
  (b) an aliphatic –CH2–OH group (R–CH2–OH) in which aside from the O atom the carbon has exactly one other carbon neighbor.
  
Improved algorithm:
  1. Parse the molecule and add explicit hydrogens.
  2. Loop over oxygen atoms to find those that look like a hydroxyl (i.e. oxygen with at least one hydrogen
     and exactly one heavy-atom neighbor). For each such “isolated –OH” we count it.
  3. For each candidate –OH, examine its single attached heavy atom (expected to be carbon). Now,
     (A) require that this carbon is sp3‐hybridized,
     (B) ignore the –OH bond and count the remaining neighbors (only heavy atoms, ignoring hydrogens) and
         count explicitly attached hydrogens.
         • If the carbon has 3 H’s and no other heavy neighbors then it is a CH3–OH candidate.
         • If the carbon has 2 H’s and exactly one heavy neighbor (i.e. an R–CH2–OH candidate) then accept it.
  4. Finally, if the molecule has a single –OH then the candidate qualifies the whole molecule.
     Otherwise, compute the ratio (candidate primary OH/total isolated OH) and require this ratio to be at least 20%.
     
  Note: This heuristic may still mis‐classify very multifunctional molecules.
"""

from rdkit import Chem

def is_primary_alcohol(smiles: str):
    """
    Determines if a molecule is classified as a primary alcohol compound based on its SMILES.
    
    A compound is considered a primary alcohol if at least one isolated -OH group is attached to a saturated carbon
    that qualifies as primary. In our method we first add explicit hydrogens, then loop over every oxygen atom to
    find –OH groups (using the criterion that the oxygen is bonded to at least one hydrogen and exactly one heavy atom).
    For each candidate –OH we examine the attached carbon’s environment: after discounting the –OH oxygen, the carbon
    should either have three hydrogens (methanol) or two hydrogens and exactly one other heavy atom (RCH2–OH).
    
    If the molecule has only one hydroxyl group and that group qualifies, it is classified as a primary alcohol.
    Otherwise, we require that at least 20% of all isolated –OH groups (those that have exactly one heavy neighbor)
    are primary.
    
    Args:
        smiles (str): SMILES string representing the molecule.
    Returns:
        bool: True if classified as a primary alcohol compound, False otherwise.
        str: Explanation of the decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens so we get reliable counts.
    mol = Chem.AddHs(mol)
    
    total_OH = 0   # count of isolated -OH groups (oxygen having >=1 H and exactly one heavy neighbor)
    primary_candidate_count = 0
    candidate_carbon_idx = None  # record index of first primary candidate carbon
    
    # Loop over atoms; we are looking for oxygen atoms representing hydroxyl groups.
    for oxygen in mol.GetAtoms():
        if oxygen.GetSymbol() != "O":
            continue
        # Count explicit hydrogens attached to the oxygen:
        h_neighbors = [nbr for nbr in oxygen.GetNeighbors() if nbr.GetSymbol() == "H"]
        if len(h_neighbors) < 1:
            continue  # no H attached
        # Collect heavy neighbors (atomic number > 1)
        heavy_neighbors = [nbr for nbr in oxygen.GetNeighbors() if nbr.GetAtomicNum() > 1]
        if len(heavy_neighbors) != 1:
            continue  # not a simple isolated hydroxyl group
        total_OH += 1
        # Get the attached heavy atom – expect it to be carbon.
        carbon = heavy_neighbors[0]
        if carbon.GetSymbol() != "C":
            continue  # not attached to carbon
        # The carbon must be sp3; this excludes cases like phenolic -OH.
        if carbon.GetHybridization() != Chem.rdchem.HybridizationType.SP3:
            continue
        
        # Now, examine the carbon environment.
        # Remove the oxygen we are currently examining from its neighbor list.
        carbon_neighbors = [nbr for nbr in carbon.GetNeighbors() if nbr.GetIdx() != oxygen.GetIdx()]
        # Count explicit hydrogens among these neighbors.
        h_count = sum(1 for nbr in carbon_neighbors if nbr.GetSymbol() == "H")
        # Count heavy neighbors (atomic number > 1)
        heavy_nb = [nbr for nbr in carbon_neighbors if nbr.GetAtomicNum() > 1]
        # For a methanol group (CH3-OH), we expect 3 hydrogens and no other heavy atom.
        if h_count == 3 and len(heavy_nb) == 0:
            primary_candidate_count += 1
            if candidate_carbon_idx is None:
                candidate_carbon_idx = carbon.GetIdx()
            continue
        # For a primary alcohol of type R-CH2-OH, we expect 2 hydrogens and exactly one heavy neighbor.
        if h_count == 2 and len(heavy_nb) == 1:
            primary_candidate_count += 1
            if candidate_carbon_idx is None:
                candidate_carbon_idx = carbon.GetIdx()
            continue
        # Otherwise, the candidate does not meet our strict primary alcohol connectivity.
    
    # If no hydroxyl groups or no candidates were found, return False.
    if total_OH == 0:
        return False, "No hydroxyl (-OH) group found"
    if primary_candidate_count == 0:
        return False, "No primary -OH (CH3-OH or RCH2-OH) group found"
    
    # Determine classification. If there is only one hydroxyl then that candidate wins.
    if total_OH == 1:
        return True, f"Contains primary alcohol group: candidate at carbon atom index {candidate_carbon_idx} (sole -OH group)."
    
    # Otherwise, use a heuristic ratio threshold.
    ratio = primary_candidate_count / total_OH
    # We lower the threshold to 20% so that compounds that have many -OH groups but only one primary one
    # (e.g. prostaglandin F2α 1-ethanolamide) are still flagged.
    threshold = 0.20
    if ratio < threshold:
        return False, f"Only {ratio*100:.1f}% of hydroxyl groups are primary; not classified as primary alcohol compound"
    
    return True, f"Contains primary alcohol group: candidate at carbon atom index {candidate_carbon_idx}. (Primary OH ratio = {ratio*100:.1f}%)"

# For testing purposes (you may remove or comment out in production)
if __name__ == '__main__':
    test_smiles = [
        "Nc1nc(=O)[nH]cc1CO",  # 5-(hydroxymethyl)cytosine (expected True)
        "CCC(CC)Nc1cc(ccc1N1C(=O)CCC1(CO)CO)C(O)=O",  # 1-[4-carboxy-2-(3-pentylamino)phenyl]-5,5-bis(hydroxymethyl)pyrrolidin-2-one (expected True)
        "[H][C@]12C[C@]3([H])C(=C)C(=O)O[C@@]3([H])[C@@H](O)\\C(C)=C\\CC\\C(C)=C\\CC[C@@]1(CO)O2",  # sinularolide C (expected True)
        "OCCC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F",  # perfluoro-1-decanol (expected True)
        "O(C([H])([H])[H])[H]",  # methanol-d1 (expected True)
        "CC(C)=CCC\\C(C)=C\\CO",  # geraniol (expected True)
        "CCCCCCCCCC=C",  # 10-undecen-1-ol (expected True)
        # Example that used to be a false negative (prostaglandin F2α 1-ethanolamide)
        "CCCCC[C@H](O)\\C=C\\[C@H]1[C@H](O)C[C@H](O)[C@@H]1C\\C=C/CCCC(=O)NCCO",
        # A couple examples from false positives (which ideally should NOT be classified as primary alcohol)
        "C[C@H]1CN(C(=O)C2=C(C=CC(=C2)NC(=O)CCCN(C)C)O[C@@H]1CN(C)C(=O)NC3CCCCC3)[C@H](C)CO",
        "C1CC(C1)CN(C[C@H]2[C@H]([C@H](N2)CO)C3=CC=CC=C3)C(=O)C4=CC=C(C=C4)F"
    ]
    for s in test_smiles:
        result, reason = is_primary_alcohol(s)
        print(f"SMILES: {s}\nResult: {result}\nReason: {reason}\n{'-'*60}")