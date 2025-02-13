"""
Classifies: CHEBI:50477 butyrate ester
"""
#!/usr/bin/env python
"""
Classifies: Butyrate ester
Definition: A butyrate ester is any carboxylic ester in which the acyl (acid) portion is derived from butyric acid.
Butyric acid (n‐butyric acid) has the structure CH3–CH2–CH2–COOH, so the acyl group becomes CH3–CH2–CH2–C(=O)–.
This program looks for an ester fragment having a linear, unbranched butyryl group.
The procedure is:
  1. Use a SMARTS pattern that explicitly captures the six heavy atoms:
     (1) a terminal CH3,
     (2) first CH2,
     (3) second CH2,
     (4) the carbonyl carbon,
     (5) the ester oxygen,
     (6) the carbonyl oxygen (double–bonded to atom 4).
  2. For each match, verify that atoms 1–3 (the butyl chain) are unbranched.
  3. Check that the carbonyl carbon (atom 4) is only attached to the chain, the ester oxygen, and the carbonyl oxygen.
  4. Verify that the carbonyl oxygen (atom 6) is only bonded to atom 4.
If these conditions hold for at least one match found, we classify the molecule as a butyrate ester.
"""

from rdkit import Chem

def is_butyrate_ester(smiles: str):
    """
    Determines if a molecule (given as a SMILES string) contains a butyrate ester group.
    A butyrate ester is defined as an ester where the acyl component is derived from 
    n-butyric acid (CH3CH2CH2C(=O)-), with a linear (unbranched) three‐carbon chain.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if a valid butyrate ester substructure is found, False otherwise.
        str: Explanation of the classification result.
    """
    # Parse the molecule from the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a SMARTS pattern that captures the butyryl (acyl) fragment.
    # Labeled as follows:
    # [CH3:1]-[CH2:2]-[CH2:3]-[C;X3:4](=[O:6])-[O;!$([O-]);!$([O+]):5]
    #   Atom1: terminal methyl (CH3)
    #   Atom2: first methylene (CH2)
    #   Atom3: second methylene (CH2)
    #   Atom4: carbonyl carbon (sp2, trigonal planar)
    #   Atom6: carbonyl oxygen (double-bonded to atom4)
    #   Atom5: ester oxygen (should be neutral)
    pattern = Chem.MolFromSmarts("[CH3:1]-[CH2:2]-[CH2:3]-[C;X3:4](=[O:6])-[O;!$([O-]);!$([O+]):5]")
    if pattern is None:
        return False, "Error in defining SMARTS pattern"

    # Get all substructure matches.
    matches = mol.GetSubstructMatches(pattern)
    if not matches:
        return False, "No butyrate ester substructure (CH3CH2CH2C(=O)O) found."

    # Function to count heavy (non-hydrogen) neighbors, excluding those in a set of allowed indices.
    def heavy_neighbor_count(atom, allowed_ids):
        return sum(1 for nb in atom.GetNeighbors() if nb.GetAtomicNum() > 1 and nb.GetIdx() not in allowed_ids)

    # Examine each match candidate.
    for match in matches:
        # Ensure we have at least six indices corresponding to our labeled atoms.
        if len(match) < 6:
            continue
        idx1, idx2, idx3, idx4, idx5, idx6 = match[:6]
        a1 = mol.GetAtomWithIdx(idx1)  # should be CH3
        a2 = mol.GetAtomWithIdx(idx2)  # CH2
        a3 = mol.GetAtomWithIdx(idx3)  # CH2
        a4 = mol.GetAtomWithIdx(idx4)  # carbonyl carbon
        a5 = mol.GetAtomWithIdx(idx5)  # ester oxygen
        a6 = mol.GetAtomWithIdx(idx6)  # carbonyl oxygen

        # Check heavy-neighbor counts to enforce unbranched linearity in the acyl chain:
        # a1 (CH3): should have exactly one heavy neighbor: a2.
        if heavy_neighbor_count(a1, {idx2}) != 0:
            continue
        # a2 (first CH2): should have exactly two heavy neighbors: a1 and a3.
        if heavy_neighbor_count(a2, {idx1, idx3}) != 0:
            continue
        # a3 (second CH2): should have exactly two heavy neighbors: a2 and a4.
        if heavy_neighbor_count(a3, {idx2, idx4}) != 0:
            continue
        # a4 (carbonyl carbon):
        # It is expected to be bonded to a3, a5 (ester oxygen) and a6 (carbonyl oxygen); no extra heavy neighbors.
        if heavy_neighbor_count(a4, {idx3, idx5, idx6}) != 0:
            continue
        # a6 (carbonyl oxygen) should be bonded only to a4.
        if heavy_neighbor_count(a6, {idx4}) != 0:
            continue
        # a5 (ester oxygen) is allowed to be bonded to additional atoms (the alkoxy part) so we do not check further.

        # If all checks pass, we have found a valid butyrate group.
        return True, ("Molecule contains a butyrate ester substructure "
                      "(CH3CH2CH2C(=O)O with a linear, unbranched acyl chain).")
    
    return False, "No unbranched butyrate ester substructure (CH3CH2CH2C(=O)O) found."

# Example usage for testing.
if __name__ == "__main__":
    test_smiles = [
        "CCCC(=O)OCC",                              # ethyl butyrate: should be True.
        "CCCC(=O)OC(C)CC",                           # sec-butyl butyrate: should be True.
        "CC(=O)OC",                                 # methyl acetate (not butyrate): False.
        "[Na+].CCCC([O-])=O",                        # sodium butyrate (acid, not ester): False.
        "O[C@@](C([N+](CC(=O)CCC)(C)C)([2H])[2H])(CC([O-])=O)[2H]",  # butyryl-L-carnitine-d3: should be True.
        "C(CCCOC(C(CC)C)=O)CC",                       # hexyl 2-methylbutanoate: might be False (branch in acyl part)
    ]
    for sm in test_smiles:
        result, reason = is_butyrate_ester(sm)
        print(f"SMILES: {sm}")
        print(f"Result: {result}\nReason: {reason}\n")