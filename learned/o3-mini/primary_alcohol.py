"""
Classifies: CHEBI:15734 primary alcohol
"""
"""
Classifies: Primary Alcohol
Definition:
  A primary alcohol is defined as a compound in which a hydroxyl (-OH) group is attached
  to a saturated carbon atom (sp3) that has either (a) three hydrogen atoms attached (CH3–OH)
  or (b) two hydrogen atoms and exactly one carbon neighbor (RCH2–OH).
  
Our improved method first adds explicit hydrogens. Then it scans every oxygen atom to see whether 
it appears to be an –OH group by checking that it has exactly one heavy-atom neighbor (typically the carbon).
For every candidate –OH group we then check the environment of its attached carbon:
  • If the carbon is sp3 and has three hydrogens then it is a methanol moiety.
  • If the carbon is sp3 and has two hydrogens and exactly one other carbon neighbor then it is an RCH2–OH.
We also count all -OH groups found. Finally, we classify the molecule as a primary alcohol if at least
one primary -OH is found and if the “primary” candidate(s) comprise at least 50% of all isolated -OH groups.
This extra filter is intended to avoid cases where an isolated primary –OH group is merely incidental within
a molecule that contains many other –OH groups in more complex environments.
    
If the molecule does not qualify or the SMILES cannot be parsed the function returns (False, reason).
"""

from rdkit import Chem

def is_primary_alcohol(smiles: str):
    """
    Determines if a molecule is a primary alcohol compound based on its SMILES string.
    
    A primary alcohol (as a compound, not just a functional group) should have at least one –OH group 
    attached to a saturated (sp3) carbon that is either CH3–OH or RCH2–OH. In addition, for our purposes
    the found primary –OH(s) should constitute a majority of the isolated –OH groups detected.
    
    Args:
        smiles (str): SMILES string representing the molecule.
        
    Returns:
        bool: True if the compound is classified as a primary alcohol, False otherwise.
        str: Explanation of the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Make sure hydrogens are explicit so that neighbor counts are reliable.
    mol = Chem.AddHs(mol)
    
    # Initialize counts and list to record candidate indices.
    primary_count = 0
    total_oh = 0
    candidate_carbon_idx = None  # to record first primary candidate carbon index
    
    # Loop over all atoms looking for potential -OH groups.
    # We require that an -OH group here is indicated by an oxygen that is bonded to at least one hydrogen and one heavy atom.
    for atom in mol.GetAtoms():
        if atom.GetSymbol() != "O":
            continue
        # Check: the oxygen should have at least one hydrogen.
        # (Here we use the explicit hydrogens which were added.)
        hs = [nbr for nbr in atom.GetNeighbors() if nbr.GetSymbol() == "H"]
        if len(hs) < 1:
            continue
        # For a simple hydroxyl group, we expect the oxygen to have exactly one heavy-atom neighbor.
        heavy_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() > 1]
        if len(heavy_neighbors) != 1:
            continue
        # Now we have a candidate –OH group.
        total_oh += 1
        carbon = heavy_neighbors[0]
        # The attached heavy atom should be carbon.
        if carbon.GetSymbol() != "C":
            continue
        # The carbon must be sp3.
        if carbon.GetHybridization() != Chem.rdchem.HybridizationType.SP3:
            continue
        
        # Count the number of hydrogens attached to the carbon.
        # (After using Chem.AddHs these counts are explicit.)
        h_count = sum(1 for nbr in carbon.GetNeighbors() if nbr.GetSymbol() == "H")
        # Also count how many carbon neighbors (excluding the oxygen we are examining)
        c_neighbors = [nbr for nbr in carbon.GetNeighbors() if nbr.GetSymbol() == "C"]
        
        # Case 1: Methanol moiety – CH3–OH (three hydrogens attached).
        if h_count == 3:
            primary_count += 1
            if candidate_carbon_idx is None:
                candidate_carbon_idx = carbon.GetIdx()
            continue
        
        # Case 2: Primary alcohol RCH2–OH: two hydrogens and exactly one carbon neighbor.
        if h_count == 2 and len(c_neighbors) == 1:
            primary_count += 1
            if candidate_carbon_idx is None:
                candidate_carbon_idx = carbon.GetIdx()
            continue
        # Else, if oxygen is bonded to a carbon that does not meet these criteria,
        # it might be a secondary or tertiary alcohol – we ignore these as candidates.
    
    # If no -OH groups were detected at all, or none of them qualified as primary:
    if total_oh == 0:
        return False, "No hydroxyl (-OH) group found"
    if primary_count == 0:
        return False, "No primary alcohol (-CH3-OH or -CH2-OH) group found"
    
    # Use a heuristic: require that at least 50% of the isolated -OH groups (with one heavy neighbor)
    # are of primary type.
    ratio = primary_count / total_oh
    if ratio < 0.5:
        return False, f"Only {ratio*100:.1f}% of hydroxyl groups are primary; not classified as primary alcohol compound"
    
    return True, f"Contains primary alcohol group: candidate at carbon atom index {candidate_carbon_idx}. (Primary OH ratio = {ratio*100:.1f}%)"


# For testing purposes (remove or comment out in production)
if __name__ == '__main__':
    test_smiles = [
        "Nc1nc(=O)[nH]cc1CO",  # 5-(hydroxymethyl)cytosine (expected True)
        "CCC(CC)Nc1cc(ccc1N1C(=O)CCC1(CO)CO)C(O)=O",  # 1-[4-carboxy-2-(3-pentylamino)phenyl]-5,5-bis(hydroxymethyl)pyrrolidin-2-one (expected True)
        "[H][C@]12C[C@]3([H])C(=C)C(=O)O[C@@]3([H])[C@@H](O)\\C(C)=C\\CC\\C(C)=C\\CC[C@@]1(CO)O2",  # sinularolide C (expected True)
        "OC1C=C(Cl)C(Cl)C=C1Cl",  # 2,4,5-trichlorocyclohexa-2,5-dien-1-ol (expected False; –OH attached to sp2 carbon or ambiguous)
        "O(C([H])([H])[H])[H]",  # methanol-d1 (expected True)
    ]
    for s in test_smiles:
        result, reason = is_primary_alcohol(s)
        print(f"SMILES: {s}\nResult: {result}\nReason: {reason}\n{'-'*40}")