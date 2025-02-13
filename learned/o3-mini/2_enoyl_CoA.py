"""
Classifies: CHEBI:19573 2-enoyl-CoA
"""
"""
Classifies: 2-enoyl-CoA
Definition: An unsaturated fatty acyl-CoA in which the S-acyl group contains a double bond between positions 2 and 3.
This program:
 - verifies the presence of a CoA moiety by looking for a characteristic CoA fragment,
 - then searches for a thioester carbonyl such that (a) it is part of a thioester (has an S neighbor)
   and (b) its acyl (alpha) carbon (the neighbor that is not S or the carbonyl oxygen) is involved in a C=C double bond.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_2_enoyl_CoA(smiles: str):
    """
    Determines if a molecule (given as a SMILES string) qualifies as a 2-enoyl-CoA.
    The molecule must have both:
      1. A CoA moiety as detected through a characteristic substructure.
      2. An acyl thioester group in which the carbonyl (position 1) is attached via S
         (making a thioester) and the next carbon (alpha carbon, position 2) is double-bonded
         to another carbon (beta carbon, position 3).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule qualifies as a 2-enoyl-CoA; otherwise False.
        str: Explanation of the classification result.
    """
    # Parse the input SMILES string into an RDKit Mol object.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."

    # First, check for the presence of a CoA moiety.
    # We use a CoA fragment that is commonly found in acyl-CoA structures.
    # This pattern represents a core portion of the pantetheine/ADP moiety.
    coa_smarts = "SCCNC(=O)CCNC(=O)"
    coa_pattern = Chem.MolFromSmarts(coa_smarts)
    if coa_pattern is None:
        return False, "Error in SMARTS for CoA fragment."
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "CoA moiety fragment not detected."

    # Now, search for a candidate thioester carbonyl whose acyl chain has a double bond between its alpha and beta carbons.
    enoyl_found = False
    # Iterate over all carbon atoms. We look for a carbonyl carbon (C=O) that is also attached to a sulfur (for thioester).
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:
            continue  # We're only interested in carbons
        neighbors = atom.GetNeighbors()
        has_carbonyl_oxygen = False
        has_thioester_sulfur = False
        # First, check bonds to determine if this carbon has a double-bonded oxygen and a single bond to S.
        for bond in atom.GetBonds():
            # Identify double bonded oxygen (carbonyl oxygen).
            if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                nbr = bond.GetOtherAtom(atom)
                if nbr.GetAtomicNum() == 8:
                    has_carbonyl_oxygen = True
            # Identify single-bonded sulfur.
            if bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                nbr = bond.GetOtherAtom(atom)
                if nbr.GetAtomicNum() == 16:
                    has_thioester_sulfur = True

        # Proceed if this atom appears to be a thioester carbonyl (has carbonyl O and an S neighbor).
        if has_carbonyl_oxygen and has_thioester_sulfur:
            # Now, among the neighbors, identify the acyl chain (the “alpha” carbon):
            # It should not be the sulfur nor the carbonyl oxygen.
            alpha_candidates = []
            for nbr in neighbors:
                if nbr.GetAtomicNum() not in (8, 16):  # not O or S
                    alpha_candidates.append(nbr)
            # Check each candidate. For a valid 2-enoyl thioester the alpha carbon (position 2)
            # must be drawn as sp2 due to the double bond with a beta carbon.
            for alpha in alpha_candidates:
                # Look for a double bond from the alpha carbon to a beta carbon.
                for bond in alpha.GetBonds():
                    # We only care about double bonds
                    if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                        beta = bond.GetOtherAtom(alpha)
                        if beta.GetAtomicNum() == 6 and beta.GetIdx() != atom.GetIdx():
                            # Found the desired alpha-beta double bond.
                            enoyl_found = True
                            break
                if enoyl_found:
                    break
        if enoyl_found:
            break

    if not enoyl_found:
        return False, "No unsaturated thioester group with a double bond between positions 2 and 3 was detected."

    return True, "Molecule contains a 2-enoyl-CoA fragment (unsaturated acyl thioester with CoA moiety)."

# Example usage (uncomment to test):
# test_smiles = "CCCCCCCCCCCCCCCCCC\\C=C\\C(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12"
# result, explanation = is_2_enoyl_CoA(test_smiles)
# print(result, explanation)