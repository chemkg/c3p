"""
Classifies: CHEBI:51689 enone
"""
"""
Classifies: enone (alpha,beta-unsaturated ketone)
An enone is defined as having the motif:
    R(1)R(2)C=CR(3)-C(=O)R(4) with R(4) ≠ H,
i.e. the carbonyl is conjugated to an alkene.
This improved definition relaxes the earlier aromatic restrictions.
"""

from rdkit import Chem

def is_enone(smiles: str):
    """
    Determines if a molecule is an enone (alpha,beta-unsaturated ketone with non-hydrogen substituent)
    based on its SMILES string.
    
    An enone is defined as a structure containing the motif:
        R(1)R(2)C=C-C(=O)-R(4)
    where R(4) is not a hydrogen. Because many valid enones have parts that are aromatic 
    (or drawn with aromatic bonds), we do not restrict the C=C atoms from being aromatic.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is an enone, False otherwise
        str: Explanation of the classification result
    """
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern for a generic enone:
    #   - "[#6]" matches any carbon atom.
    #   - "= [#6]" matches a double bond between two carbons.
    #   - "-[#6](=O)" matches a carbon single bonded to the alkene carbon and doubly bonded to an oxygen.
    #   - "-[!#1]" ensures that the substituent bonded to the carbonyl carbon is not hydrogen.
    # This gives a motif: C=C-C(=O)-R, where R is any heavy atom.
    enone_smarts = "[#6]=[#6]-[#6](=O)-[!#1]"
    enone_pattern = Chem.MolFromSmarts(enone_smarts)
    if enone_pattern is None:
        return False, "Error in SMARTS pattern"
    
    # Find all substructure matches for the enone motif.
    matches = mol.GetSubstructMatches(enone_pattern)
    if not matches:
        return False, "No enone motif (C=C-C(=O)-R with R ≠ H) found in the molecule"
    
    # Verify candidate matches by ensuring that the carbonyl carbon is indeed double-bonded to an oxygen.
    # Each match is a 4-tuple: (alpha, beta, carbonyl, substituent)
    for match in matches:
        if len(match) != 4:
            # Skip if the match is not exactly 4 atoms.
            continue
        alpha_idx, beta_idx, carbonyl_idx, substituent_idx = match
        
        # Get carbonyl atom.
        carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
        
        # Check if the carbonyl atom is double bonded to at least one oxygen.
        found_double_bonded_oxygen = False
        for neighbor in carbonyl_atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 8:  # oxygen
                bond = mol.GetBondBetweenAtoms(carbonyl_atom.GetIdx(), neighbor.GetIdx())
                if bond is not None and bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                    found_double_bonded_oxygen = True
                    break
        if not found_double_bonded_oxygen:
            # Candidate rejected – the ketone group is not valid.
            continue
        
        # If we reach here, we found a valid enone motif.
        return True, "Enone motif found: C=C-C(=O)-R (with R ≠ H)"
    
    # If no candidate passes filtering, then no valid enone motif is found.
    return False, "No valid enone motif found after filtering candidate matches"

# Example usage and testing:
if __name__ == "__main__":
    test_smiles = [
        "O(C=1C(=C(O)C=C(OC)C1)C(=O)/C=C/C2=CC=C(OC)C=C2)",  # Flavokawain A, expected True
        "[H][C@@]12CCC(=O)[C@@]1(C)CC[C@@]1([H])[C@@]2([H])[C@H](O)CC2=CC(=O)CC[C@]12C",  # Steroid enone, expected True
        "O(C1=C(CC=C(C)C)C(O)=C(C=C1)C(=O)/C=C/C2=CC=C(O)C=C2)C",  # 4-Hydroxyderricin, expected True
        "O=C(/C=C/C)C",  # 3E-Nonen-2,5-dione, expected True
        "O=C(/C=C/CC)C",  # 6-Methyl-3-hepten-2-one, trans-, expected True
    ]
    for s in test_smiles:
        result, explanation = is_enone(s)
        print(f"SMILES: {s}\nResult: {result}\nExplanation: {explanation}\n")