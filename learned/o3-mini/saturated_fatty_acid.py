"""
Classifies: CHEBI:26607 saturated fatty acid
"""
"""
Classifies: Saturated fatty acid

A saturated fatty acid is defined here as a molecule that:
  - Contains only C, H, and O atoms.
  - Contains exactly one carboxylic acid group defined as [CX3](=O)[OX2H].
  - Has no carbon–carbon multiple bonds.
  - Contains no additional carbonyl groups (C=O bonds between C and O) outside the acid group if the molecule has 7 or more carbons.
  - Does not have an -OH substituent on the alpha carbon (the carbon directly bonded to the acid carbon).
  - May have a single additional hydroxyl group at the end of the chain (the ω‐position). Any extra oxygen-containing functionality (e.g. non‐terminal –OH, keto, methoxy) is disallowed.
  
Note: Unlike the previous attempt, we do not reject molecules solely on the basis of ring structures.
"""

from rdkit import Chem

def is_saturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a saturated fatty acid based on its SMILES string.
    The molecule must:
      - Contain only C, H, and O atoms.
      - Have exactly one carboxylic acid group (SMARTS: [CX3](=O)[OX2H]).
      - Have no carbon–carbon double or triple bonds.
      - In molecules with 7 or more carbons, have no additional carbonyl groups (C=O) beyond that in the acid.
      - Not have an -OH substituent on the alpha carbon (the carbon directly bonded to the acid carbon).
      - Optionally have one extra (terminal) hydroxyl group at the opposite omega (ω) end.
      
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule is classified as a saturated fatty acid, False otherwise.
        str: Explanation/reason for the decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # (1) Allowed atoms: only H (1), C (6), and O (8)
    allowed_atomic_nums = {1, 6, 8}
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in allowed_atomic_nums:
            return False, "Contains atoms other than C, H, and O"
    
    # (2) Identify the carboxylic acid group.
    acid_smarts = "[CX3](=O)[OX2H]"
    acid_pattern = Chem.MolFromSmarts(acid_smarts)
    acid_matches = mol.GetSubstructMatches(acid_pattern)
    if not acid_matches:
        return False, "No carboxylic acid group found"
    if len(acid_matches) != 1:
        return False, "More than one carboxylic acid group found"
    
    # Record acid group atom indices.
    acid_match = acid_matches[0]
    acid_carbon_idx = acid_match[0]       # the carboxyl carbon
    acid_carbonyl_oxygen_idx = acid_match[1]  # the allowed carbonyl oxygen
    acid_oh_oxygen_idx = acid_match[2]     # the hydroxyl oxygen of the COOH
    
    # (3) Check for any carbon-to-carbon multiple bonds.
    for bond in mol.GetBonds():
        if bond.GetBondType() in [Chem.rdchem.BondType.DOUBLE, Chem.rdchem.BondType.TRIPLE]:
            a1 = bond.GetBeginAtom()
            a2 = bond.GetEndAtom()
            if a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 6:
                return False, "Contains carbon-to-carbon multiple bonds"
    
    # (4) Count carbons in the molecule.
    total_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    
    # (5) Check for additional carbonyl groups (C=O bonds between C and O, not part of the acid group).
    # For molecules with 7 or more carbons, any extra C=O (double bonds: C=O) will lead to rejection.
    if total_carbons >= 7:
        for bond in mol.GetBonds():
            if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                a1 = bond.GetBeginAtom()
                a2 = bond.GetEndAtom()
                # We look for a C=O bond:
                if (a1.GetAtomicNum(), a2.GetAtomicNum()) in [(6,8), (8,6)]:
                    idxs = {bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()}
                    # Allow if this is the acid carbonyl
                    if {acid_carbon_idx, acid_carbonyl_oxygen_idx} == idxs:
                        continue
                    else:
                        return False, "Contains additional carbonyl group(s) beyond the carboxylic acid"
    
    # (6) Check that the acid group is terminal in the sense that there is no -OH on the alpha carbon.
    acid_carbon = mol.GetAtomWithIdx(acid_carbon_idx)
    for neighbor in acid_carbon.GetNeighbors():
        # Skip the neighbor that is the carbonyl oxygen.
        if neighbor.GetIdx() == acid_carbonyl_oxygen_idx:
            continue
        # For the remaining neighbor (alpha carbon), check if it bears an -OH group.
        for subnbr in neighbor.GetNeighbors():
            if subnbr.GetAtomicNum() == 8:
                bond = mol.GetBondBetweenAtoms(neighbor.GetIdx(), subnbr.GetIdx())
                if bond and bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                    # If the oxygen has at least one explicit hydrogen, treat it as an -OH.
                    if subnbr.GetTotalNumHs() > 0:
                        return False, "Contains an alpha-hydroxyl group adjacent to the carboxylic acid"
    
    # (7) Check for any extra oxygen functionality (outside the acid group).
    # Allowed extra oxygen: one terminal hydroxyl group (at the ω-end), defined as an -OH on a carbon that is terminal (has exactly one neighboring carbon).
    # Gather oxygen indices that are not part of the acid group.
    acid_atom_idxs = {acid_carbon_idx, acid_carbonyl_oxygen_idx, acid_oh_oxygen_idx}
    extra_oxygen_idxs = [atom.GetIdx() for atom in mol.GetAtoms() 
                           if atom.GetAtomicNum() == 8 and atom.GetIdx() not in acid_atom_idxs]
    
    allowed_terminal_oh_count = 0
    for oidx in extra_oxygen_idxs:
        oxy = mol.GetAtomWithIdx(oidx)
        # We expect these oxygens to be in a hydroxyl group (i.e. bonded by a single bond to a carbon and having at least one H)
        if oxy.GetDegree() != 1 or oxy.GetTotalNumHs() < 1:
            return False, "Contains extra oxygen functionality that is not a terminal hydroxyl group"
        # Get the carbon to which this oxygen is attached.
        nbr = oxy.GetNeighbors()[0]
        if nbr.GetAtomicNum() != 6:
            return False, "Extra oxygen is not attached to a carbon"
        # Check that the carbon is terminal (only one carbon neighbor among heavy atoms).
        carbon_neighbors = [a for a in nbr.GetNeighbors() if a.GetAtomicNum() == 6]
        if len(carbon_neighbors) != 1:
            return False, "Contains an extra hydroxyl group on a non-terminal carbon"
        allowed_terminal_oh_count += 1
        if allowed_terminal_oh_count > 1:
            return False, "More than one terminal hydroxyl group found"
    
    # Passed all checks.
    return True, ("Molecule has exactly one carboxylic acid group, only C/H/O atoms, no carbon–carbon multiple bonds, "
                  "no disallowed additional carbonyl or alpha‐OH groups, and at most one terminal extra hydroxyl group")

# Example usage (you may add more SMILES strings to test):
if __name__ == '__main__':
    test_smiles = [
        "OCCCCCCCCCCCCCCCCCCCC(O)=O",      # 20-hydroxyicosanoic acid (should be True)
        "CC(C)C(C)C(O)=O",                 # 2,3-dimethylbutyric acid (should be True)
        "CC(=O)CCC(O)=O",                  # 4-oxopentanoic acid (allowed for a short chain)
        "CCCCCCCCCCCCCCCC(O)=O",           # tridecanoic acid (should be True)
        "C[C@@H](O)CCCCCCCCC[C@@H](O)CC(O)=O",  # Example that should be rejected (extra –OH on non-terminal carbons)
        "CCCCCCCC[C@H]1CCC[C@@H]1CCCCCCC(O)=O",  # lactobacillic acid variant with ring (should be True now)
    ]
    for sm in test_smiles:
        result, reason = is_saturated_fatty_acid(sm)
        print(f"SMILES: {sm}\n -> {result}, Reason: {reason}\n")