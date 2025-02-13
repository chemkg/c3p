"""
Classifies: CHEBI:16749 1-phosphatidyl-1D-myo-inositol
"""
"""
Classifies: 1-phosphatidyl-1D-myo-inositol
Definition: A phosphatidylinositol in which the inositol moiety is the 1D-myo isomer
            (i.e. a cyclohexane carrying six oxygens with five free hydroxyl groups, with the 
            phosphate substituting the remaining hydroxyl) and the phosphatidyl (diacylglycerol) 
            unit provides at least two acyl (ester) groups.
The function is_1_phosphatidyl_1D_myo_inositol takes a SMILES string and returns a boolean and 
a reason for the classification.
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_1_phosphatidyl_1D_myo_inositol(smiles: str):
    """
    Determines whether a molecule (given as a SMILES string) is a 1-phosphatidyl-1D-myo-inositol.
    
    The new strategy is broken down into two parts:
    
    1. Check for the inositol-phosphate motif:
       - Add explicit hydrogens so we can verify hydroxyl groups.
       - Look for a phosphate center (P(=O)(O)(O)) and then consider each oxygen neighbor
         attached to that phosphorus.
       - For each such oxygen, see if its non-phosphate neighbor is part of a 6-membered ring 
         whose atoms are all carbons.
       - For that ring, count how many ring carbons have a free hydroxyl group (an O–H substituent).
         In 1D-myo-inositol with a phosphate at position 1 the ring should provide five “OH” groups.
    
    2. Check for the diacylglycerol unit:
       - In phosphatidylinositols the glycerol part is acylated on at least two positions.
       - We use an ester SMARTS "C(=O)O[C]" to find ester bonds representing acyl chains.
       - To avoid catching the phosphate-ester, we discard any ester in which the bridging O is 
         also bound to a phosphorus.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as 1-phosphatidyl-1D-myo-inositol, False otherwise.
        str: A reason detailing the decision.
    """
    
    # Parse SMILES, and add explicit hydrogens (to help count OH groups)
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    mol = Chem.AddHs(mol)
    
    # --- Part 1: Look for the inositol-phosphate moiety ---
    # Define a SMARTS for a phosphate group: P(=O)(O)(O)
    phos_smarts = "P(=O)(O)(O)"
    phos_query = Chem.MolFromSmarts(phos_smarts)
    if phos_query is None:
        return False, "Error in phosphate SMARTS pattern."
    
    phos_matches = mol.GetSubstructMatches(phos_query)
    valid_inositol_found = False
    inositol_reason = ""
    
    # Get ring information (list of rings: each ring is a tuple of atom indices)
    ring_info = mol.GetRingInfo().AtomRings()
    
    # For each P match, examine its neighbors
    for match in phos_matches:
        # The match gives the index of the phosphorus atom.
        phos_idx = match[0]
        phos_atom = mol.GetAtomWithIdx(phos_idx)
        # Iterate over neighbors of phosphorus; we want a neighbor that is oxygen but not another P.
        for neigh in phos_atom.GetNeighbors():
            if neigh.GetAtomicNum() != 8:
                continue
            # Skip if this oxygen is also connected to another phosphorus (could be part of a polyphosphate)
            if any(n.GetAtomicNum() == 15 for n in neigh.GetNeighbors() if n.GetIdx() != phos_idx):
                continue
            # For this oxygen, find the other neighbor (besides the phosphorus)
            nonP_neighbors = [n for n in neigh.GetNeighbors() if n.GetIdx() != phos_idx]
            if not nonP_neighbors:
                continue
            donor_atom = nonP_neighbors[0]
            # We now want donor_atom to be part of a 6-membered ring composed solely of carbon.
            candidate_rings = []
            for ring in ring_info:
                if donor_atom.GetIdx() in ring and len(ring) == 6:
                    # Check that all atoms in this ring are carbons.
                    if all(mol.GetAtomWithIdx(idx).GetAtomicNum() == 6 for idx in ring):
                        candidate_rings.append(ring)
            if not candidate_rings:
                continue
            # For each candidate ring, count how many ring carbons bear an -OH (oxygen with at least one H)
            for ring in candidate_rings:
                oh_count = 0
                for idx in ring:
                    carbon = mol.GetAtomWithIdx(idx)
                    # Look for neighbors of the ring carbon that are not in the ring.
                    for nb in carbon.GetNeighbors():
                        if nb.GetIdx() in ring:
                            continue
                        # Check if neighbor is oxygen and has at least one hydrogen attached.
                        if nb.GetAtomicNum() == 8 and nb.GetTotalNumHs() >= 1:
                            oh_count += 1
                            break  # count one –OH per carbon max.
                # In a 1D-myo configuration with phosphate replacement at one position, expect 5 OH groups.
                if oh_count == 5:
                    valid_inositol_found = True
                    break
            if valid_inositol_found:
                break
        if valid_inositol_found:
            break
    
    if not valid_inositol_found:
        return False, "No valid inositol moiety with phosphate at position 1 (6-membered carbon ring with five -OH groups) was found."
    
    # --- Part 2: Look for diacylglycerol acyl (ester) units ---
    # Use a simple SMARTS for an ester bond: a carbonyl (C(=O)) attached to an O which is attached to a carbon.
    # In the SMARTS, atom indices: 0: carbonyl C, 1: carbonyl O, 2: bridging O, 3: next carbon.
    ester_smarts = "C(=O)O[C]"
    ester_query = Chem.MolFromSmarts(ester_smarts)
    if ester_query is None:
        return False, "Error in ester SMARTS pattern."
    
    ester_matches = mol.GetSubstructMatches(ester_query)
    acyl_count = 0
    for match in ester_matches:
        # The bridging oxygen is match index 2.
        bridging_O = mol.GetAtomWithIdx(match[2])
        # Exclude this ester if the bridging oxygen is attached to any phosphorus (i.e. belonging to phosphate linkage)
        if any(neigh.GetAtomicNum() == 15 for neigh in bridging_O.GetNeighbors()):
            continue
        acyl_count += 1
    
    if acyl_count < 2:
        return False, f"Found {acyl_count} acyl ester group(s); at least 2 are expected for a diacylglycerol unit."
    
    return True, "Molecule matches 1-phosphatidyl-1D-myo-inositol structural criteria."

# (Optional) Testing the function on one provided SMILES:
if __name__ == "__main__":
    # Example test: one of the provided valid SMILES
    test_smiles = "[C@](COC(=O)CCCCCCC/C=C\\C/C=C\\CCCC)(COP(O[C@H]1C(C(C([C@H](C1O)O)O)O)O)(O)=O)(OC(=O)CCCCC/C=C\\C/C=C\\C/C=C\\C/C=C\\CCCCC)[H]"
    result, reason = is_1_phosphatidyl_1D_myo_inositol(test_smiles)
    print("Result:", result)
    print("Reason:", reason)