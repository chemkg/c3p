"""
Classifies: CHEBI:50477 butyrate ester
"""
"""
Classifies: Butyrate ester
Definition: Any carboxylic ester for which the carboxylic acid component is butyric acid.
Butyric acid (n‐butyric acid) has structure CH3CH2CH2COOH, so in the ester its acyl group becomes CH3CH2CH2C(=O)–.
We look for an ester fragment defined as CH3–CH2–CH2–C(=O)–O (with the ester oxygen neutral)
and we post–check that the three-carbon chain (CH3, CH2, CH2) is unbranched.
"""

from rdkit import Chem

def is_butyrate_ester(smiles: str):
    """
    Determines if a molecule (given as a SMILES string) is a butyrate ester.
    A butyrate ester is defined as a carboxylic ester having the acyl group of n-butyric acid,
    namely CH3CH2CH2C(=O)–, and the chain must be linear (unbranched), with a neutral ester oxygen.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if a butyrate ester substructure (CH3CH2CH2C(=O)O with an unbranched acyl chain) is found,
              False otherwise.
        str: Explanation of the classification result.
    """
    # Parse the molecule from the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern for an ester fragment where the acyl component is butyric acid:
    # We label the atoms as:
    # [#6:1] is the terminal methyl (CH3),
    # [#6:2] and [#6:3] are the methylene groups (CH2),
    # [C:4](=O) is the carbonyl carbon, 
    # [O;!$([O-]);!$([O+]):5] is the ester oxygen (the oxygen must not be charged).
    #
    # The pattern therefore is: CH3-CH2-CH2-C(=O)-O, with labels 1–5.
    pat = Chem.MolFromSmarts("[#6:1]-[#6:2]-[#6:3]-[C:4](=O)-[O;!$([O-]);!$([O+]):5]")
    if pat is None:
        return False, "Error in defining SMARTS pattern"
    
    # Get all substructure matches.
    matches = mol.GetSubstructMatches(pat)
    if not matches:
        return False, "No butyrate ester substructure (CH3CH2CH2C(=O)O) found."
    
    # Function to check if an atom in the chain is unbranched.
    # For our chain we expect:
    # - The terminal CH3 (atom1) to only have one heavy-atom neighbor (atom2).
    # - The first CH2 (atom2) to have exactly two heavy-atom neighbors (atoms 1 and 3).
    # - The second CH2 (atom3) to have exactly two heavy-atom neighbors (atoms 2 and the carbonyl, atom4).
    def is_unbranched(atom, allowed_neighbor_ids):
        # Check heavy atoms (atomic number > 1).
        for nb in atom.GetNeighbors():
            if nb.GetAtomicNum() > 1 and nb.GetIdx() not in allowed_neighbor_ids:
                return False
        return True

    # Check each match candidate from the SMARTS.
    for match in matches:
        # Ensure we only take the first five indices corresponding to our labeled atoms.
        # (Sometimes the match tuple returns extra indices.)
        if len(match) < 5:
            continue
        idx1, idx2, idx3, idx4, idx5 = match[:5]
        atom1 = mol.GetAtomWithIdx(idx1)
        atom2 = mol.GetAtomWithIdx(idx2)
        atom3 = mol.GetAtomWithIdx(idx3)
        # For atom1 (CH3), allow only atom2.
        if not is_unbranched(atom1, {idx2}):
            continue
        # For atom2 (first CH2), allow only atoms 1 and 3.
        if not is_unbranched(atom2, {idx1, idx3}):
            continue
        # For atom3 (second CH2), allow only atoms 2 and 4.
        if not is_unbranched(atom3, {idx2, idx4}):
            continue
        # If all extra checks pass, we have identified an unbranched butyrate ester group.
        return True, "Molecule contains a butyrate ester substructure (CH3CH2CH2C(=O)O) with a linear acyl chain."
    
    # If no candidate match passes the unbranched test, return false.
    return False, "No unbranate (linear) butyrate ester substructure (CH3CH2CH2C(=O)O) found."

# Example usage for testing.
if __name__ == "__main__":
    # A few test SMILES strings:
    test_smiles = [
        "CCCC(=O)OCC",         # ethyl butyrate: should return True.
        "CCCC(=O)OC(C)CC",      # sec-butyl butyrate: should return True.
        "CC(=O)OC",            # methyl acetate (not butyrate): should return False.
        "[Na+].CCCC([O-])=O",   # sodium butyrate (acid, not ester): should return False.
        "O[C@@](C([N+](CC(=O)CCC)(C)C)([2H])[2H])(CC([O-])=O)[2H]",  # butyryl-L-carnitine-d3: should return True.
    ]
    for sm in test_smiles:
        result, reason = is_butyrate_ester(sm)
        print(f"SMILES: {sm}")
        print(f"Result: {result}, Reason: {reason}\n")