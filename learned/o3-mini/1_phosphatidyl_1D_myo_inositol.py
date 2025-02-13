"""
Classifies: CHEBI:16749 1-phosphatidyl-1D-myo-inositol
"""
"""
Classifies: 1-phosphatidyl-1D-myo-inositol
Definition: A phosphatidylinositol in which the inositol moiety is the 1D-myo isomer 
 (i.e. a cyclohexane ring whose five carbons display free –OH groups while the remaining position 
 carries a phospho substituent), and the phosphatidyl (diacylglycerol) unit provides at least two acyl (ester) groups.
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_1_phosphatidyl_1D_myo_inositol(smiles: str):
    """
    Determines whether a molecule (given as a SMILES string) is a 1-phosphatidyl-1D-myo-inositol.
    
    This function performs two major checks:
    
    1. Inositol moiety check:
       - The molecule is first parsed and explicit hydrogens are added.
       - Every ring of size 6 that is composed solely of carbons (i.e. cyclohexane cores) is examined.
       - For each such ring, for every carbon we inspect atoms attached (that are not in the ring).
         If an oxygen is directly attached:
             • If that oxygen in turn is bound to a phosphorus, we count that as the phospho substituent.
             • Otherwise, if the oxygen carries at least one hydrogen, we count that as a free -OH.
       - A valid 1D-myo-inositol is detected if exactly one carbon shows a phospho substituent (indicating the phosphate at position 1)
         and the other five carbons show a free hydroxyl group.
    
    2. Diacylglycerol ester unit check:
       - Using a SMARTS pattern for esters, we look for acyl groups whose bridging oxygen is NOT attached to any phosphorus.
       - At least 2 acyl ester groups are required.
       
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if classified as 1-phosphatidyl-1D-myo-inositol; otherwise False.
        str: Reason detailing the decision.
    """
    # Parse SMILES and add explicit hydrogens to correctly see hydroxyl groups.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    mol = Chem.AddHs(mol)
    
    # --- Part 1: Validate the inositol ring with a phosphate at one position ---
    ring_info = mol.GetRingInfo().AtomRings()
    valid_inositol_found = False
    inositol_reason = ""
    
    # Iterate over all rings and look for 6-membered rings
    for ring in ring_info:
        if len(ring) != 6:
            continue
        # Check that all atoms in the ring are carbons.
        if not all(mol.GetAtomWithIdx(idx).GetAtomicNum() == 6 for idx in ring):
            continue
        
        # For this cyclohexane, check substituents.
        phospho_count = 0  # count of substituents that are O attached to P
        oh_count = 0       # count of substituents that are free hydroxyl (O with H)
        
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            # Find neighbors not in the ring
            substituents = [nbr for nbr in atom.GetNeighbors() if nbr.GetIdx() not in ring]
            # For each carbon, we allow a maximum of one relevant substituent.
            found_sub = False
            for sub in substituents:
                if sub.GetAtomicNum() != 8:
                    continue
                # Check if this oxygen is attached to a phosphorus (phospho substituent)
                attached_P = any(nbr.GetAtomicNum() == 15 for nbr in sub.GetNeighbors() if nbr.GetIdx() != idx)
                if attached_P:
                    phospho_count += 1
                    found_sub = True
                    break  # only count one per carbon
                else:
                    # Check if oxygen has at least one hydrogen (free OH)
                    # GetTotalNumHs gives the count of explicit H attached.
                    if sub.GetTotalNumHs() >= 1:
                        oh_count += 1
                        found_sub = True
                        break
            # Optionally, if no substituent found, we ignore.
        # For a valid 1D-myo-inositol moiety, expect exactly one phospho substituent and five free –OH groups.
        if phospho_count == 1 and oh_count == 5:
            valid_inositol_found = True
            break  # stop after finding one valid inositol ring

    if not valid_inositol_found:
        return False, "No valid inositol moiety with phosphate at position 1 (6-membered cyclohexane with 5 free –OH groups and 1 phospho substituent) was found."
    
    # --- Part 2: Validate the diacylglycerol (acyl ester) units ---
    # SMARTS for an ester bond: carbonyl (C(=O)) bound to an oxygen that in turn is bound to a carbon.
    ester_smarts = "C(=O)O[C]"
    ester_query = Chem.MolFromSmarts(ester_smarts)
    if ester_query is None:
        return False, "Error in ester SMARTS pattern."
    
    ester_matches = mol.GetSubstructMatches(ester_query)
    acyl_count = 0
    for match in ester_matches:
        # In the match the bridging oxygen is at index 2.
        bridging_O = mol.GetAtomWithIdx(match[2])
        # Exclude ester bonds if the bridging oxygen is attached to a phosphorus (i.e. part of the phosphate link).
        if any(neigh.GetAtomicNum() == 15 for neigh in bridging_O.GetNeighbors()):
            continue
        acyl_count += 1
    
    if acyl_count < 2:
        return False, f"Found {acyl_count} acyl ester group(s); at least 2 are expected for the diacylglycerol unit."
    
    return True, "Molecule matches 1-phosphatidyl-1D-myo-inositol structural criteria."

# (Optional) Testing the function on one example valid SMILES:
if __name__ == "__main__":
    # One of the known valid SMILES examples:
    test_smiles = "[C@](COC(=O)CCCCCCC/C=C\\C/C=C\\CCCC)(COP(O[C@H]1C(C(C([C@H](C1O)O)O)O)O)(O)=O)(OC(=O)CCCCC/C=C\\C/C=C\\C/C=C\\C/C=C\\CCCCC)[H]"
    result, reason = is_1_phosphatidyl_1D_myo_inositol(test_smiles)
    print("Result:", result)
    print("Reason:", reason)