"""
Classifies: CHEBI:16749 1-phosphatidyl-1D-myo-inositol
"""
"""
Classifies: 1-phosphatidyl-1D-myo-inositol
Definition: A phosphatidylinositol in which the inositol moiety is the 1D-myo isomer 
 (i.e. a cyclohexane ring in which five carbons display a free –OH group while the remaining 
 carbon is substituted with an oxygen that is bound to a phosphate group) and the phosphatidyl unit 
 (diacylglycerol) provides at least two acyl (ester) groups.
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_1_phosphatidyl_1D_myo_inositol(smiles: str):
    """
    Determines whether a molecule (given as a SMILES string) is a 1-phosphatidyl-1D-myo-inositol.
    
    The function performs two major checks:
    
    1. Inositol moiety check:
       - After parsing and adding explicit hydrogens, each six-membered ring that is composed 
         entirely of carbon atoms is examined. For a valid 1D-myo inositol, every carbon in 
         the ring must have exactly one oxygen neighbor (a substituent). Among these 6 oxygens, 
         exactly one must be bonded to a phosphorus atom (representing the phosphate substituent) 
         and the other five must be free hydroxyl groups (i.e. each oxygen has at least one attached hydrogen).
         
    2. Diacylglycerol (ester) unit check:
       - Using a SMARTS pattern for ester bonds (C(=O)O[C]), count the acyl groups that are not 
         part of a phosphate linkage. At least two acyl ester groups must be present.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as 1-phosphatidyl-1D-myo-inositol, False otherwise.
        str: Explanation for the classification.
    """
    # Parse SMILES string and add explicit hydrogens
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    mol = Chem.AddHs(mol)
    
    # --- Part 1: Look for a valid inositol ring ---
    ring_info = mol.GetRingInfo().AtomRings()
    valid_inositol_found = False
    
    # Iterate over all rings of size 6
    for ring in ring_info:
        if len(ring) != 6:
            continue
        # Check that every atom in this ring is a carbon (atomic number 6)
        if not all(mol.GetAtomWithIdx(idx).GetAtomicNum() == 6 for idx in ring):
            continue
        
        # For each ring carbon, we expect exactly one substituent oxygen
        phospho_count = 0  # count of carbons with oxygen attached to phosphorus (phospho substituent)
        oh_count = 0       # count of carbons with a free –OH (oxygen with at least one hydrogen)
        valid_ring = True  # assume valid until proven otherwise
        
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            # Get neighbors not in the ring
            substituents = [nbr for nbr in atom.GetNeighbors() if nbr.GetIdx() not in ring]
            # We expect exactly one oxygen substituent on each ring carbon.
            oxygen_subs = [nbr for nbr in substituents if nbr.GetAtomicNum() == 8]
            if len(oxygen_subs) != 1:
                valid_ring = False
                break
            # Check the lone oxygen substituent
            o_atom = oxygen_subs[0]
            # Determine if this oxygen is bound to phosphorus (i.e. phosphate substituent)
            attached_to_P = any(neigh.GetAtomicNum() == 15 for neigh in o_atom.GetNeighbors() if neigh.GetIdx() != idx)
            if attached_to_P:
                phospho_count += 1
            else:
                # For free hydroxyl, check that the oxygen has at least one hydrogen attached.
                if o_atom.GetTotalNumHs() >= 1:
                    oh_count += 1
                else:
                    valid_ring = False
                    break
        # For a valid 1D-myo inositol, we expect exactly 6 oxygen substituents,
        # exactly one phospho substituent and five free –OH groups.
        if valid_ring and phospho_count == 1 and oh_count == 5:
            valid_inositol_found = True
            break

    if not valid_inositol_found:
        return False, ("No valid inositol moiety with a phosphate substituent at one position was found. "
                       "Expected a cyclohexane ring with six oxygen substituents (one bound to phosphorus and "
                       "five as free –OH groups).")
    
    # --- Part 2: Check for the diacylglycerol (ester) unit ---
    # Define a SMARTS pattern for an ester bond: a carbonyl (C(=O)) linked to an oxygen that is in turn linked to a carbon.
    ester_smarts = "C(=O)O[C]"
    ester_query = Chem.MolFromSmarts(ester_smarts)
    if ester_query is None:
        return False, "Error in ester SMARTS pattern."
    
    ester_matches = mol.GetSubstructMatches(ester_query)
    acyl_count = 0
    for match in ester_matches:
        # match indices: e.g. match[0] = carbonyl carbon, match[1] = carbonyl oxygen, match[2] = bridging carbon
        bridging_O = mol.GetAtomWithIdx(match[1])
        # Exclude this ester if the bridging oxygen is attached to phosphorus (phosphate bond)
        if any(neigh.GetAtomicNum() == 15 for neigh in bridging_O.GetNeighbors()):
            continue
        acyl_count += 1
        
    if acyl_count < 2:
        return False, f"Found {acyl_count} acyl ester group(s); at least 2 are required for the diacylglycerol moiety."
    
    return True, "Molecule matches the 1-phosphatidyl-1D-myo-inositol structural criteria."

# (Optional) Testing the function on a known example:
if __name__ == "__main__":
    test_smiles = "[C@](COC(=O)CCCCCCC/C=C\\C/C=C\\CCCC)(COP(O[C@H]1C(C(C([C@H](C1O)O)O)O)O)(O)=O)(OC(=O)CCCCC/C=C\\C/C=C\\C/C=C\\C/C=C\\CCCCC)[H]"
    result, reason = is_1_phosphatidyl_1D_myo_inositol(test_smiles)
    print("Result:", result)
    print("Reason:", reason)