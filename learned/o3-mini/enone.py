"""
Classifies: CHEBI:51689 enone
"""
"""
Classifies: an enone (alpha,beta-unsaturated ketone with R(4) ≠ H)
An enone must have a conjugated C=C–C(=O) moiety.
For example, an enone follows the general formula:
    R(1)R(2)C=CR(3)-C(=O)R(4)
with the requirement that R(4) is not hydrogen.
In this implementation, we require that the two carbons forming the C=C 
are non-aromatic (to avoid mis‐matching tightly aromatic or quinone-like systems)
and that the carbonyl carbon is bound to a heavy atom (ensuring R(4) ≠ H).
"""

from rdkit import Chem

def is_enone(smiles: str):
    """
    Determines if a molecule is an enone based on its SMILES string.
    
    An enone is defined as an alpha,beta-unsaturated ketone where the carbonyl is in conjugation
    with a carbon–carbon double bond and has a substituent (R(4) ≠ H) at the carbonyl carbon.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if the molecule is an enone, False otherwise
        str: Explanation/reason for the classification
    """
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Revised SMARTS explanation:
    # We look for a four-atom fragment corresponding to: 
    #   C(alpha)=C(beta)-C(=O)-R4
    # where we require the following:
    # 1) The two carbons in the double bond are non-aromatic: [C;!a]
    # 2) The carbonyl carbon (the third atom) must be a trigonal carbon attached to a carbonyl oxygen.
    # 3) The fourth atom (R4) is any heavy atom (implicitly ensuring R4 is not hydrogen).
    #
    # This SMARTS may not catch every edge case but tends to avoid matches
    # in complex fused aromatic systems.
    enone_smarts = "[C;!a]=[C;!a]-[C](=O)[#6]"
    enone_pattern = Chem.MolFromSmarts(enone_smarts)
    if enone_pattern is None:
        return False, "Error in SMARTS pattern"
    
    # Find all matches for the enone fragment
    matches = mol.GetSubstructMatches(enone_pattern)
    if not matches:
        return False, "No enone motif (conjugated C=C–C(=O) with R(4) ≠ H) found in the molecule"
    
    # Post-filter each candidate match:
    # For each four-atom fragment match, we can check:
    # - That the carbonyl carbon (third atom) indeed has a double-bonded oxygen neighbor.
    # - That the double bond is conjugated by verifying that the bond between the beta carbon and 
    #   the carbonyl carbon is a double bond.
    # Note: Some false positives appear when the motif occurs as part of an aromatic fused system.
    # Here we already restricted the double bond atoms to be non-aromatic.
    for match in matches:
        a_idx, b_idx, c_idx, d_idx = match
        a = mol.GetAtomWithIdx(a_idx)  # alpha carbon
        b = mol.GetAtomWithIdx(b_idx)  # beta carbon
        c = mol.GetAtomWithIdx(c_idx)  # carbonyl carbon
        d = mol.GetAtomWithIdx(d_idx)  # substituent attached to carbonyl (R4)
        
        # Verify that the carbonyl carbon actually is doubly bonded to an O
        has_carbonyl_oxygen = False
        for nb in c.GetNeighbors():
            if nb.GetAtomicNum() == 8:  # oxygen
                bond = mol.GetBondBetweenAtoms(c.GetIdx(), nb.GetIdx())
                if bond is not None and bond.GetBondTypeAsDouble() == 2.0:
                    has_carbonyl_oxygen = True
                    break
        if not has_carbonyl_oxygen:
            continue  # not a valid ketone motif

        # Check that the bond between b and c is a double bond.
        bond_bc = mol.GetBondBetweenAtoms(b.GetIdx(), c.GetIdx())
        if bond_bc is None or bond_bc.GetBondTypeAsDouble() != 2.0:
            continue  # the conjugation is missing

        # If we reached here, we have found an enone fragment.
        # (Additional filtering could be added here if needed.)
        return True, "Enone motif found: alpha,beta-unsaturated ketone with non-hydrogen substituent on the carbonyl"
    
    # If no candidate passed the post-filters, return False.
    return False, "No valid enone motif found after filtering candidate matches"

# Example usage:
if __name__ == "__main__":
    test_smiles = [
        "O(C=1C(=C(O)C=C(OC)C1)C(=O)/C=C/C2=CC=C(OC)C=C2)C",  # Flavokawain A, expected True
        "[H][C@@]12CCC(=O)[C@@]1(C)CC[C@@]1([H])[C@@]2([H])[C@H](O)CC2=CC(=O)CC[C@]12C",  # a steroid enone, expected True
        "C[C@H]1O[C@H](C[C@@H](O)[C@@H]1O)c1ccc2C(=O)C3=C([C@H](O)C[C@]4(O)C[C@@](C)(O)CC(=O)[C@]34O)C(=O)c2c1O",  # False positive example from previous attempt
        "OC1=C([C@]2([C@@H](CCC(=C2)C)C(C)C)[H])C(O)=C(C(O)=C1[C@]3([C@@H](CCC(=C3)C)C(C)C)[H])C(=O)CCC4=CC=CC=C4"  # (-)-Neolinderatin, expected True
    ]
    for s in test_smiles:
        res, reason = is_enone(s)
        print(f"SMILES: {s}\nResult: {res}\nReason: {reason}\n")