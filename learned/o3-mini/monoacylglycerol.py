"""
Classifies: CHEBI:17408 monoacylglycerol
"""
#!/usr/bin/env python
"""
Classifies: monoacylglycerol
Definition:
  A glyceride in which any one of the R groups (position not specified) is an acyl group 
  while the remaining two R groups can be either H or alkyl groups.

Heuristic:
  1. Parse the SMILES string.
  2. Look for a three‐carbon (propane) chain with the expected pattern for glycerol:
     the first and third carbon (CH2) and the central carbon (CH).
  3. For each such chain, examine the substituents (neighbors not in the chain). 
     A monoacylglycerol should have exactly three oxygen substituents – one ester linkage 
     (–O–C(=O)–) and two free hydroxyls (–OH).
  4. If found, return True with a reason; otherwise return False with an explanation.
"""

from rdkit import Chem

def is_monoacylglycerol(smiles: str):
    """
    Determines if a molecule is a monoacylglycerol based on its SMILES string.
    A monoacylglycerol is a glyceride in which one of the three substituents on a glycerol
    backbone (CH2-CHOH-CH2) is acylated (via an ester bond) and the other two groups are free
    (as hydroxyls or simple alkyl, though in our heuristic we look for -OH groups).
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is classified as a monoacylglycerol, False otherwise.
        str: Explanation/reason for the classification.
    """
    # Attempt to parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # --- Helper functions ---
    def is_ester_oxygen(oxy):
        """
        Return True if the oxygen atom 'oxy' is part of an ester bond.
        We assume that if an oxygen bonded to a carbon that itself has a double bond to another oxygen,
        then that oxygen is esterified.
        """
        # Iterate over neighbors of the oxygen atom.
        for nbr in oxy.GetNeighbors():
            # We expect the neighbor to be carbon.
            if nbr.GetAtomicNum() == 6:
                # Look at bonds from the carbon.
                for bond in nbr.GetBonds():
                    # Skip the bond that connects back to oxy.
                    if bond.GetOtherAtom(nbr) == oxy:
                        continue
                    # If the bond is a double bond and the other atom is oxygen, we likely have C(=O)
                    if bond.GetBondTypeAsDouble() == 2 and bond.GetOtherAtom(nbr).GetAtomicNum() == 8:
                        return True
        return False

    def is_free_hydroxyl(oxy):
        """
        Return True if the oxygen atom 'oxy' appears to be a free hydroxyl.
        We expect a hydroxyl oxygen to have at least one hydrogen attached.
        """
        # Many times the explicit H count is not present; use GetTotalNumHs.
        if oxy.GetTotalNumHs() > 0:
            return True
        return False

    # --- Step 1: Look for a glycerol backbone candidate
    # We look for a three-carbon chain (CH2-CH-CH2).
    # Terminal carbons are CH2 (degree=2 for heteroatoms when ignoring the chain) and central is CH.
    glycerol_smarts = "[CH2;D2]-[CH;D3]-[CH2;D2]"
    glycerol_pattern = Chem.MolFromSmarts(glycerol_smarts)
    matches = mol.GetSubstructMatches(glycerol_pattern)
    if not matches:
        return False, "No three-carbon (glycerol) backbone found"

    # Examine each candidate backbone for the expected oxygen substituents.
    for match in matches:
        # match is a tuple of atom indices for the 3 carbons in order.
        # For a glycerol (or monoacylglycerol) backbone, each carbon should have one extra oxygen
        # substituent (apart from the carbon–carbon bonds connecting the three carbons).
        ester_count = 0
        free_oh_count = 0
        ok_candidate = True

        # Loop over the three backbone carbons.
        for idx in match:
            atom = mol.GetAtomWithIdx(idx)
            # Find neighbors that are not in the backbone (i.e. substituents).
            subst_neigh = [nbr for nbr in atom.GetNeighbors() if nbr.GetIdx() not in match]
            # We expect exactly one substituent oxygen for each glycerol carbon.
            # (Note: In monoacylglycerol, one of these oxygens will be esterified.)
            oxy_neighbors = [nbr for nbr in subst_neigh if nbr.GetAtomicNum() == 8]
            if len(oxy_neighbors) != 1:
                ok_candidate = False
                break
            # Check the oxygen substituent.
            oxy = oxy_neighbors[0]
            if is_ester_oxygen(oxy):
                ester_count += 1
            elif is_free_hydroxyl(oxy):
                free_oh_count += 1
            else:
                # It is an oxygen but not clearly -OH, so mark candidate as not matching.
                ok_candidate = False
                break

        # For monoacylglycerol we expect one ester and two free hydroxyl groups.
        if ok_candidate:
            if ester_count == 1 and free_oh_count == 2:
                return True, "Glycerol backbone found with one acyl (ester) substituent and two free hydroxyl groups"
            else:
                # Candidate found but its substituent pattern does not match monoacylglycerol.
                # Continue to check any other candidate.
                continue

    return False, "No glycerol backbone with exactly one acyl (ester) substitution and two free hydroxyls found"

    
# Example usage:
if __name__ == '__main__':
    # Test a few examples; you can add more test SMILES.
    test_smiles = [
        "CCCCCCCCCCCC(=O)OC[C@@H](O)CO",  # 1-hexadecanoyl-sn-glycerol (monoacylglycerol)
        "CCCC(=O)OCC(O)CO",              # monobutyrin (also monoacylglycerol)
        "OCC(O)CO"                      # glycerol (no ester group)
    ]
    for s in test_smiles:
        res, reason = is_monoacylglycerol(s)
        print(f"SMILES: {s}\n Result: {res}\n Reason: {reason}\n")