"""
Classifies: CHEBI:87658 decanoate ester
"""
"""
Classifies: Decanoate ester

A decanoate ester is defined as the fatty acid ester derived from decanoic (capric) acid.
Decanoic acid has the structure CH3(CH2)8C(=O)OH so that its ester (decanoate) acyl fragment is
CH3(CH2)8C(=O)O (or its deprotonated form). This classifier looks for that decanoate fragment
(with exactly 10 carbons in the chain, including the carbonyl carbon) and then verifies that:
  - the terminal CH3 is not further substituted (degree = 1),
  - the ester oxygen bridges to another group (degree >= 2),
to ensure we are not matching free decanoic acid.
"""
from rdkit import Chem

def is_decanoate_ester(smiles: str):
    """
    Determines if a molecule contains a decanoate ester moiety (formed via condensation of 
    decanoic acid with an alcohol or phenol) based on its SMILES string.
    
    The decanoate moiety built from decanoic acid (CH3(CH2)8C(=O)OH) in the ester form is:
      CH3-(CH2)8-C(=O)-O  (or with a negative charge on the oxygen)
    
    To capture this, we search for a substructure defined by the two SMARTS patterns:
      Pattern 1 (protonated): "[CH3][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2]C(=O)O"
      Pattern 2 (deprotonated): "[CH3][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2]C(=O)[O-]"
    
    After a match is found, we verify:
       - The terminal CH3 atom (first atom in match) has degree exactly 1.
       - The ester oxygen (last atom) is bridging (degree >=2), to ensure it links to another R group.
       
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if a valid decanoate ester fragment is found, False otherwise.
        str: Explanation for the classification.
    """
    # Parse the SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the SMARTS patterns for the decanoate acyl fragment.
    # Note: We do not include the ;D1 attribute here and will check the degree explicitly.
    pattern_prot = "[CH3][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2]C(=O)O"
    pattern_deprot = "[CH3][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2]C(=O)[O-]"
    
    smarts_patterns = [pattern_prot, pattern_deprot]
    
    for pat in smarts_patterns:
        substruct = Chem.MolFromSmarts(pat)
        if substruct is None:
            continue  # Should not occur, but in case the SMARTS is invalid
        
        matches = mol.GetSubstructMatches(substruct)
        if not matches:
            continue  # Try the next pattern if no match found
        
        # Process each match; we expect exactly 11 atoms in the matched substructure:
        # Index 0: CH3
        # Indices 1-8: eight CH2 groups
        # Index 9: carbonyl carbon (C=O)
        # Index 10: ester oxygen
        for match in matches:
            if len(match) != 11:
                continue  # Not an exact match for our decanoate fragment
            
            # Verify the terminal CH3 is truly terminal: its degree must equal 1.
            ch3_atom = mol.GetAtomWithIdx(match[0])
            if ch3_atom.GetDegree() != 1:
                # The purported terminal methyl has additional connections.
                continue
            
            # Verify the ester oxygen is bridging to an external group:
            oxy_atom = mol.GetAtomWithIdx(match[10])
            if oxy_atom.GetDegree() < 2:
                # Likely matching free decanoic acid (or not attached to a substituent)
                continue
            
            # Found a valid decanoate ester fragment.
            return True, "Contains decanoate ester moiety with a terminal decanoate acyl chain"
    
    # If no match passes the tests, return False with an explanation.
    return False, "Decanoate ester fragment not found or does not appear as a terminal decanoate acyl chain"

# Example usage (for testing purposes):
if __name__ == "__main__":
    test_smiles = [
        "[C@H](OC(CCCCCCCCC)=O)(C[N+](C([2H])([2H])[2H])(C)C)CC(=O)[O-]",
        "O=C(OC)CCCCCCCCC",
        "CCCCCCCCCC(=O)OCC(COP(O)(O)=O)OC(=O)CCCCCCCCC"
    ]
    for s in test_smiles:
        result, reason = is_decanoate_ester(s)
        print(f"SMILES: {s}\nResult: {result}\nReason: {reason}\n")