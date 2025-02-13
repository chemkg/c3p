"""
Classifies: CHEBI:51689 enone
"""
"""
Classifies: enone (alpha,beta-unsaturated ketone, where the carbonyl substituent is not hydrogen)
An enone is defined as having the motif:
    R(1)R(2)C=C-C(=O)R(4)  (with R(4) ≠ H)
Here we require that the two atoms in the C=C double bond and the carbonyl carbon are not aromatic.
The carbonyl carbon must have a double-bonded oxygen (verified by the SMARTS) and the substituent (R(4))
must be a heavy atom (i.e. not hydrogen).
"""

from rdkit import Chem

def is_enone(smiles: str):
    """
    Determines if a molecule is an enone (alpha,beta-unsaturated ketone with non-hydrogen substituent)
    based on its SMILES string.

    An enone is defined as a structure containing C(alpha)=C(beta)-C(=O)-R,
    where the two carbons involved in the C=C bond and the carbonyl carbon are non-aromatic,
    the carbonyl carbon is doubly bonded to an oxygen, and R (R(4)) is a non-hydrogen atom.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an enone, False otherwise
        str: Explanation for the classification result
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a SMARTS pattern for the enone motif.
    # [C;!a]=[C;!a] captures a non-aromatic C=C double bond.
    # [C;!a](=O) ensures the next carbon (the carbonyl carbon) is non-aromatic and has a double bond to oxygen.
    # -[!#1] specifies that the substituent on the carbonyl is not hydrogen.
    enone_smarts = "[C;!a]=[C;!a]-[C;!a](=O)-[!#1]"
    enone_pattern = Chem.MolFromSmarts(enone_smarts)
    if enone_pattern is None:
        return False, "Error in SMARTS pattern"

    # Get all substructure matches for the enone motif.
    matches = mol.GetSubstructMatches(enone_pattern)
    if not matches:
        return False, "No enone motif (C=C-C(=O)-R with R ≠ H) found in the molecule"
    
    # Iterate through candidate matches (each should be a 4-tuple)
    for match in matches:
        if len(match) != 4:
            # Skip if not exactly 4 atoms are matched.
            continue
        # Unpack the match: positions of alpha, beta, carbonyl, and substituent atoms
        alpha_idx, beta_idx, carbonyl_idx, substituent_idx = match
        alpha_atom = mol.GetAtomWithIdx(alpha_idx)
        beta_atom = mol.GetAtomWithIdx(beta_idx)
        carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
        substituent_atom = mol.GetAtomWithIdx(substituent_idx)
        
        # (Optional extra filtering:)
        # Verify that the carbonyl atom in the match indeed is doubly bonded to an oxygen.
        has_double_bonded_oxygen = False
        for nb in carbonyl_atom.GetNeighbors():
            if nb.GetAtomicNum() == 8:  # oxygen
                bond = mol.GetBondBetweenAtoms(carbonyl_atom.GetIdx(), nb.GetIdx())
                # Check if the bond is double
                if bond is not None and bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                    has_double_bonded_oxygen = True
                    break
        if not has_double_bonded_oxygen:
            continue  # This candidate does not have the proper carbonyl group
        
        # If reached here, we expect a valid enone motif.
        return True, "Enone motif found: C(alpha)=C(beta)-C(=O)-R (with non-hydrogen substituent)"
    
    # If no match passed extra filtering, return False.
    return False, "No valid enone motif found after filtering candidate matches"

# Example usage and testing:
if __name__ == "__main__":
    test_smiles = [
        "O(C=1C(=C(O)C=C(OC)C1)C(=O)/C=C/C2=CC=C(OC)C=C2)",  # Flavokawain A, expected True
        "[H][C@@]12CCC(=O)[C@@]1(C)CC[C@@]1([H])[C@@]2([H])[C@H](O)CC2=CC(=O)CC[C@]12C",  # Steroid enone, expected True
        "C[C@H]1O[C@H](C[C@@H](O)[C@@H]1O)c1ccc2C(=O)C3=C([C@H](O)C[C@]4(O)C[C@@](C)(O)CC(=O)[C@]34O)C(=O)c2c1O",  # False positive candidate previously
        "OC1=C([C@]2([C@@H](CCC(=C2)C)C(C)C)[H])C(O)=C(C(O)=C1[C@]3([C@@H](CCC(=C3)C)C(C)C)[H])C(=O)CCC4=CC=CC=C4"  # (-)-Neolinderatin, expected True
    ]
    for s in test_smiles:
        result, explanation = is_enone(s)
        print(f"SMILES: {s}\nResult: {result}\nExplanation: {explanation}\n")