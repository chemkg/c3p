"""
Classifies: CHEBI:50477 butyrate ester
"""
"""
Classifies: Butyrate ester (Any carboxylic ester where the carboxylic acid component is butyric acid)

Butyric acid is n-butyric acid (CH3CH2CH2COOH). In an ester the acyl group becomes CH3CH2CH2C(=O)-.
We search for a fragment corresponding to a linear, unbranched n-butyrate moiety attached
via an ester bond. To address pitfalls of the previous approach, we define a SMARTS for an ester fragment
with a 3-carbon chain preceding the carbonyl. Then we post-check the candidate match to ensure that the chain atoms 
are unbranched (i.e. are only connected to each other via the appropriate bonds) and that the ester oxygen is neutral.
"""

from rdkit import Chem

def is_butyrate_ester(smiles: str):
    """
    Determines if a molecule (given as a SMILES string) is a butyrate ester.
    A butyrate ester is defined as a carboxylic ester having the acyl group of n-butyric acid,
    that is CH3CH2CH2C(=O)-. To help avoid false positives we require that the three-carbon chain is linear,
    and that the ester oxygen is not charged.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if a butyrate ester substructure (i.e. CH3CH2CH2C(=O)O with a linear acyl chain) is found,
              False otherwise.
        str: Explanation of the classification result.
    """
    # Parse the SMILES into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a SMARTS pattern for an ester fragment with a four-carbon acyl group.
    # We label the atoms as follows:
    # [#6:1]-[#6:2]-[#6:3]-[C:4](=O)-[O;!$([O-]);!$([O+]):5]
    # Here atom1 should be the terminal methyl (CH3), atom2 and atom3 are methylene groups,
    # atom4 is the carbonyl carbon, and atom5 is the connected oxygen (which must be neutral).
    pat = Chem.MolFromSmarts("[#6:1]-[#6:2]-[#6:3]-[C:4](=O)-[O;!$([O-]);!$([O+]):5]")
    if pat is None:
        return False, "Error in defining SMARTS pattern"
    
    matches = mol.GetSubstructMatches(pat)
    if not matches:
        return False, "No butyrate ester substructure (CH3CH2CH2C(=O)O) found."

    # Function to check if a chain atom is unbranched.
    # For the chain we expect:
    # - Atom1 (terminal methyl) should have only one heavy-atom neighbor (namely atom2)
    # - Atom2 (the first methylene) should have exactly 2 heavy-atom neighbors (atom1 and atom3)
    # - Atom3 (the second methylene) should have exactly 2 heavy-atom neighbors (atom2 and atom4)
    # Here heavy atoms are those with atomic number > 1.
    def is_unbranched(atom, allowed_neighbor_ids):
        for nb in atom.GetNeighbors():
            if nb.GetAtomicNum() > 1 and nb.GetIdx() not in allowed_neighbor_ids:
                return False
        return True

    # Now check each match candidate to ensure the chain is unbranched.
    for match in matches:
        # match is a tuple of atom indices corresponding to [atom1, atom2, atom3, atom4, atom5]
        idx1, idx2, idx3, idx4, idx5 = match
        atom1 = mol.GetAtomWithIdx(idx1)
        atom2 = mol.GetAtomWithIdx(idx2)
        atom3 = mol.GetAtomWithIdx(idx3)
        # For atom1, allowed neighbor is only atom2.
        if not is_unbranched(atom1, {idx2}):
            continue
        # For atom2, allowed neighbors are atom1 and atom3.
        if not is_unbranched(atom2, {idx1, idx3}):
            continue
        # For atom3, allowed neighbors are atom2 and atom4.
        if not is_unbranched(atom3, {idx2, idx4}):
            continue
        # If we have passed the branching test, we have found an unbranched n-butyrate acyl group.
        return True, "Molecule contains an ester group with the n-butyric acid moiety (CH3CH2CH2C(=O)O) with a linear acyl chain."

    # If no match candidate passed the extra filtering, none qualifies as an n-butyrate ester.
    return False, "No unbranched butyrate ester substructure (CH3CH2CH2C(=O)O) found."

# Example usage (testing some SMILES)
if __name__ == "__main__":
    test_smiles = [
        "CCCC(=O)OCC",         # ethyl butyrate: should return True.
        "CCCC(=O)OC(C)CC",      # sec-butyl butyrate: should return True.
        "CC(=O)OC",            # methyl acetate (not butyrate): should return False.
        "[Na+].CCCC([O-])=O",   # sodium butyrate (not an ester): should return False.
        "O[C@@](C([N+](CC(=O)CCC)(C)C)([2H])[2H])(CC([O-])=O)[2H]",  # butyryl-L-carnitine-d3 (should return True)
    ]
    for sm in test_smiles:
        result, reason = is_butyrate_ester(sm)
        print(f"SMILES: {sm}")
        print(f"Result: {result}, Reason: {reason}\n")