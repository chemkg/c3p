"""
Classifies: CHEBI:47622 acetate ester
"""
"""
Classifies: Any carboxylic ester where the carboxylic acid component is acetic acid.
Here we decide a molecule contains an acetate ester if it includes a bond sequence
O–C(=O)–CH3. That is, an oxygen atom bonded (via a single bond) to a carbon which is
double bonded to an oxygen and whose other substituent is a methyl group.
"""

from rdkit import Chem

def is_acetate_ester(smiles: str):
    """
    Determines if a molecule given by a SMILES string contains an acetate ester group.
    An acetate ester is defined as a carboxylic ester where the acyl part is acetic acid,
    i.e. it contains the motif –O–C(=O)–CH3.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule contains an acetate ester moiety, False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse SMILES into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # We iterate over all bonds looking for a single bond between an oxygen and a carbon.
    # In an ester bond, the oxygen (ester O) is single-bonded to a carbon that is also
    # double-bonded to an oxygen (the C=O). Furthermore, the other substituent on the carbonyl
    # should be a methyl group (i.e. only attached to the carbonyl carbon).
    for bond in mol.GetBonds():
        # We look only at single bonds.
        if bond.GetBondType() != Chem.rdchem.BondType.SINGLE:
            continue
            
        # Get the atoms at both ends of this bond.
        atom1 = bond.GetBeginAtom()
        atom2 = bond.GetEndAtom()
        
        # We require one atom to be oxygen and the other to be carbon.
        if (atom1.GetAtomicNum() == 8 and atom2.GetAtomicNum() == 6) or (atom2.GetAtomicNum() == 8 and atom1.GetAtomicNum() == 6):
            # Identify which is O and which is C.
            if atom1.GetAtomicNum() == 8:
                oxy = atom1
                c = atom2
            else:
                oxy = atom2
                c = atom1
            
            # Check that the carbon has a double bonded oxygen (i.e. it is a carbonyl carbon)
            carbonyl_found = False
            for nb in c.GetNeighbors():
                # Skip the oxygen already connected in the ester bond.
                if nb.GetIdx() == oxy.GetIdx():
                    continue
                # Check if the bond is double and if the neighbor is oxygen.
                bond_to_nb = c.GetBondBetweenAtoms(c.GetIdx(), nb.GetIdx())
                if bond_to_nb is not None and bond_to_nb.GetBondType() == Chem.rdchem.BondType.DOUBLE and nb.GetAtomicNum() == 8:
                    carbonyl_found = True
                    carbonyl_oxygen = nb
                    break
            if not carbonyl_found:
                continue  # Not an ester carbonyl, go on
            
            # Now, for an ester, the carbonyl carbon (c) must have two substituents:
            # one is the oxygen (already considered) and the other is the acyl substituent.
            # We want that substituent to be a methyl group.
            # Find the other neighbor (other than the oxygen from the O–C bond).
            acyl_candidate = None
            for nb in c.GetNeighbors():
                if nb.GetIdx() in [oxy.GetIdx(), carbonyl_oxygen.GetIdx()]:
                    continue
                acyl_candidate = nb
                break
            if acyl_candidate is None:
                continue
            
            # Check that the candidate is a carbon.
            if acyl_candidate.GetAtomicNum() != 6:
                continue
            
            # In a methyl group, the only heavy-atom connection should be to the carbonyl.
            # (Note: implicit hydrogens are not counted in GetDegree().)
            if acyl_candidate.GetDegree() != 1:
                continue
            
            # If we get here, we have found a bond O–C with c being carbonyl and its other neighbor being CH3.
            return True, "Contains an acetate ester moiety (O-C(=O)-CH3)"
            
    return False, "No acetate ester substructure (O-C(=O)-CH3) found"
    
# Example usage:
if __name__ == "__main__":
    # Some tests; the provided examples can be used to validate classification.
    test_smiles_list = [
        # A few examples from the provided list:
        ("CC(=O)Oc1ccc(cc1)C(O)=O", "4-acetoxy benzoic acid (should be acetate ester)"),
        ("CCOC(=O)C(C)C#N", "ethyl 2-cyanopropionate (belongs to class?)"),
    ]
    
    for smi, name in test_smiles_list:
        result, reason = is_acetate_ester(smi)
        print(f"SMILES: {smi}\n{name} => {result}: {reason}\n")