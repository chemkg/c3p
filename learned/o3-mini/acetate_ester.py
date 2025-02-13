"""
Classifies: CHEBI:47622 acetate ester
"""
"""
Classifies: Any carboxylic ester where the carboxylic acid component is acetic acid.
We define an acetate ester as a molecule containing the motif -O-C(=O)-CH3.
This program uses RDKit for pattern matching and substructure checking.
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
    # Parse the SMILES into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Iterate over all bonds to search for ester bonds with the acetate motif.
    for bond in mol.GetBonds():
        # Consider only single bonds
        if bond.GetBondType() != Chem.rdchem.BondType.SINGLE:
            continue

        # Get atoms at both ends of the bond
        atom1 = bond.GetBeginAtom()
        atom2 = bond.GetEndAtom()
        
        # One atom must be oxygen and the other must be carbon.
        if (atom1.GetAtomicNum() == 8 and atom2.GetAtomicNum() == 6) or (atom2.GetAtomicNum() == 8 and atom1.GetAtomicNum() == 6):
            # Identify which is the oxygen and which is the carbon
            if atom1.GetAtomicNum() == 8:
                oxy = atom1
                c = atom2
            else:
                oxy = atom2
                c = atom1
            
            # Check that the carbon (c) has a double-bonded oxygen (i.e., a carbonyl group)
            carbonyl_found = False
            carbonyl_oxygen = None
            for nb in c.GetNeighbors():
                # Skip the oxygen already connected in the ester bond.
                if nb.GetIdx() == oxy.GetIdx():
                    continue
                # Use the Mol object to get the bond between the carbon and the neighbor
                bond_to_nb = mol.GetBondBetweenAtoms(c.GetIdx(), nb.GetIdx())
                if bond_to_nb is not None and bond_to_nb.GetBondType() == Chem.rdchem.BondType.DOUBLE and nb.GetAtomicNum() == 8:
                    carbonyl_found = True
                    carbonyl_oxygen = nb
                    break
            if not carbonyl_found:
                continue  # No carbonyl found on this bond; skip
            
            # Now find the substituent on c that is not the ester oxygen (oxy) or the carbonyl oxygen.
            # For an acetate ester, this substituent (the acyl part) should be a methyl group.
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
            
            # In a methyl group, the only heavy-atom connection should be to the carbonyl carbon.
            # (Implicit hydrogens are not counted in GetDegree(), so the degree should be 1.)
            if acyl_candidate.GetDegree() != 1:
                continue
            
            # All conditions are satisfied: we have an O–C single bond where the C is double-bonded to an O 
            # and the extra substituent is a methyl group. Thus, the structure contains an acetate ester.
            return True, "Contains an acetate ester moiety (O-C(=O)-CH3)"
    
    # If none of the bonds match the criteria, return False.
    return False, "No acetate ester substructure (O-C(=O)-CH3) found"

# Example usage with test cases (for development illustration):
if __name__ == "__main__":
    test_smiles_list = [
        ("CC(=O)Oc1ccc(cc1)C(O)=O", "4-acetoxy benzoic acid (should be acetate ester)"),
        ("CCOC(=O)C(C)C#N", "ethyl 2-cyanopropionate (unlikely an acetate ester)"),
        ("CCCCCCCCOC(C)=O", "nonyl acetate (should be acetate ester)"),
    ]
    
    for smi, description in test_smiles_list:
        result, reason = is_acetate_ester(smi)
        print(f"SMILES: {smi}\n{description} => {result}: {reason}\n")