"""
Classifies: CHEBI:16749 1-phosphatidyl-1D-myo-inositol
"""
"""
Classifies: 1-phosphatidyl-1D-myo-inositol
Definition: A phosphatidylinositol in which the inositol moiety is the 1D-myo isomer
            and the phosphatidyl group is located at its position 1.
The function is_1_phosphatidyl_1D_myo_inositol takes a SMILES string and returns a boolean
and a reason for the classification.
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_1_phosphatidyl_1D_myo_inositol(smiles: str):
    """
    Determines whether a molecule (given as a SMILES string) is a 1-phosphatidyl-1D-myo-inositol.
    
    The approach here is to check three structural features:
    1. The molecule must contain an inositol moiety having the arrangement of hydroxyl groups
       characteristic of 1D-myo-inositol. (We use a SMARTS pattern that matches a cyclohexane
       ring bearing six hydroxyl groups with stereochemistry typical for D-myo-inositol.)
    2. One of the hydroxyl groups (the one corresponding to position 1) must be substituted
       (acylated) by a phosphate group. In other words, an oxygen of the inositol ring must be
       directly bonded to a phosphorus atom having the typical P(=O)(O)O fragment.
    3. In phosphatidylinositols the phosphate group is connected to a diacylglycerol unit;
       here we search for at least two ester (acyl) groups (i.e. a C(=O)O– fragment not directly
       linked to P) as a proxy for the acyl chains.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is classified as 1-phosphatidyl-1D-myo-inositol, False otherwise.
        str: A reason message explaining the classification result.
    """
    
    # Parse the SMILES string to a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # Define a SMARTS pattern for 1D-myo-inositol.
    # This pattern covers a cyclohexane ring with a hydroxyl group on each vertex.
    # The stereochemistry tags in the pattern help target the 1D-myo isomer.
    inositol_smarts = "O[C@@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O"
    inositol_mol = Chem.MolFromSmarts(inositol_smarts)
    if inositol_mol is None:
        return False, "Error in inositol SMARTS pattern."
    
    inositol_matches = mol.GetSubstructMatches(inositol_mol)
    if not inositol_matches:
        return False, "No inositol moiety with the expected 1D-myo configuration was found."
    
    # Check that at least one of the inositol matches has an oxygen that is connected to a phosphorus.
    phosphate_attached = False
    for match in inositol_matches:
        # In our SMARTS pattern the first atom (index 0) is the oxygen substituent that (in PI)
        # should be connected to the phosphate.
        inositol_O_idx = match[0]
        atom_O = mol.GetAtomWithIdx(inositol_O_idx)
        # Look at neighbors of this oxygen to see if one is a phosphorus (atomic number 15)
        for neighbor in atom_O.GetNeighbors():
            if neighbor.GetAtomicNum() == 15:
                phosphate_attached = True
                break
        if phosphate_attached:
            break
                
    if not phosphate_attached:
        return False, "Inositol moiety is not substituted with a phosphate group at position 1."
    
    # Next, look for ester bonds that represent acyl chains on the diacylglycerol part.
    # We use a simple SMARTS for an ester bond: a carbonyl (C=O) connected to an oxygen
    # that is then bound to a carbon.
    ester_smarts = "C(=O)O[C]"
    ester_mol = Chem.MolFromSmarts(ester_smarts)
    if ester_mol is None:
        return False, "Error in ester SMARTS pattern."
    ester_matches = mol.GetSubstructMatches(ester_mol)
    
    # Filter out ester matches that belong to the phosphate linkage.
    # We do so by demanding that the oxygen in the ester (atom with index corresponding to 'O' in the pattern)
    # is not also a neighbor of a phosphorus.
    acyl_count = 0
    for match in ester_matches:
        # In our pattern, atom indices: 0 is carbonyl C, 1 is carbonyl O, 2 is the bridging O, 3 is next C.
        bridging_O_idx = match[2]
        bridging_O = mol.GetAtomWithIdx(bridging_O_idx)
        bonded_P = any(neigh.GetAtomicNum() == 15 for neigh in bridging_O.GetNeighbors())
        if not bonded_P:
            acyl_count += 1
    
    if acyl_count < 2:
        return False, f"Found {acyl_count} acyl ester group(s); at least 2 are expected for a diacylglycerol unit."
    
    # If all checks pass, we classify the molecule as a 1-phosphatidyl-1D-myo-inositol.
    return True, "Molecule matches 1-phosphatidyl-1D-myo-inositol structural criteria."

# (Optional) Below lines could be used for testing.
if __name__ == "__main__":
    # Example test (one of the provided valid SMILES)
    test_smiles = "[C@](COC(=O)CCCCCCC/C=C\\C/C=C\\CCCC)(COP(O[C@H]1C(C(C([C@H](C1O)O)O)O)O)(O)=O)(OC(=O)CCCCC/C=C\\C/C=C\\C/C=C\\C/C=C\\CCCCC)[H]"
    result, reason = is_1_phosphatidyl_1D_myo_inositol(test_smiles)
    print("Result:", result)
    print("Reason:", reason)