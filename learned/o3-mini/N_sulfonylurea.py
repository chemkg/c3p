"""
Classifies: CHEBI:76983 N-sulfonylurea
"""
"""
Classifies: N-sulfonylurea
A urea in which one of the hydrogens attached to a nitrogen of the urea group is replaced by a sulfonyl group.
This moiety is common in various herbicides and antidiabetic drugs.
"""

from rdkit import Chem

def is_N_sulfonylurea(smiles: str):
    """
    Determines if a molecule is a N-sulfonylurea based on its SMILES string.
    A N-sulfonylurea is defined as a urea (N-C(=O)-N) in which one of the hydrogens on a nitrogen is replaced
    by a sulfonyl group (-S(=O)(=O)-R). In some cases an intervening amino group may be present linking the urea
    nitrogen to the sulfonyl moiety.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is N-sulfonylurea, False otherwise
        str: Reason for the classification
    """
    # Parse the input SMILES.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Use a SMARTS pattern to find a urea core.
    # This pattern looks for an sp3 (or slightly unsaturated) nitrogen, a carbonyl carbon, and another nitrogen.
    urea_smarts = "[NX3]C(=O)[NX3]"
    urea_core = Chem.MolFromSmarts(urea_smarts)
    if urea_core is None:
        return False, "Could not interpret urea SMARTS pattern"
    
    urea_matches = mol.GetSubstructMatches(urea_core)
    if not urea_matches:
        return False, "No urea (N-C(=O)-N) core found in the molecule"
    
    def is_sulfonyl(sulfur_atom):
        """
        Checks whether a sulfur atom looks like a sulfonyl center.
        We require that it has exactly two oxygen neighbors connected via a double bond and that its total degree is 4.
        """
        dbl_oxygen_count = 0
        for bond in sulfur_atom.GetBonds():
            # Get the neighbor atom.
            neighbor = bond.GetOtherAtom(sulfur_atom)
            if neighbor.GetAtomicNum() == 8 and bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                dbl_oxygen_count += 1
        # For a typical sulfonyl S, count of doubly bound oxygens is exactly 2
        # and the total degree (number of connected atoms) is 4.
        if dbl_oxygen_count == 2 and sulfur_atom.GetDegree() == 4:
            return True
        return False

    def sulfonyl_on_neighbor(n_atom, core_atoms, visited=None):
        """
        Checks whether the given nitrogen atom (n_atom) has a substituent that is (or leads to) a sulfonyl S.
        Only considers bonds that are not part of the urea core (whose indices are given in core_atoms).
        Allows one intervening nitrogen.
        visited: set of atom indices already visited to avoid loops.
        """
        if visited is None:
            visited = set()
        visited.add(n_atom.GetIdx())
        for nbr in n_atom.GetNeighbors():
            # Skip if this neighbor is part of the urea core.
            if nbr.GetIdx() in core_atoms:
                continue
            # Direct attachment: if the neighbor is sulfur and bonded by a single bond.
            if nbr.GetAtomicNum() == 16:
                bond = mol.GetBondBetweenAtoms(n_atom.GetIdx(), nbr.GetIdx())
                if bond is not None and bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                    if is_sulfonyl(nbr):
                        return True
            # Allow one intervening amino group.
            if nbr.GetAtomicNum() == 7 and nbr.GetIdx() not in visited:
                # Check that the bond from n_atom to this intervening nitrogen is single.
                bond = mol.GetBondBetweenAtoms(n_atom.GetIdx(), nbr.GetIdx())
                if bond is None or bond.GetBondType() != Chem.rdchem.BondType.SINGLE:
                    continue
                # Now, from the intervening nitrogen, check its neighbors (except going back to n_atom
                # and anything from the urea core).
                for nbr2 in nbr.GetNeighbors():
                    if nbr2.GetIdx() == n_atom.GetIdx() or nbr2.GetIdx() in core_atoms:
                        continue
                    if nbr2.GetAtomicNum() == 16:
                        bond2 = mol.GetBondBetweenAtoms(nbr.GetIdx(), nbr2.GetIdx())
                        if bond2 is not None and bond2.GetBondType() == Chem.rdchem.BondType.SINGLE:
                            if is_sulfonyl(nbr2):
                                return True
        return False

    # For every urea core match we found, check the two terminal nitrogens.
    # In our urea SMARTS "N-C(=O)-N", the match tuple has indices:
    # match[0] -> first N, match[1] -> carbonyl C, match[2] -> second N.
    for match in urea_matches:
        core_atoms = set(match)  # the three atoms forming the urea core
        n1 = mol.GetAtomWithIdx(match[0])
        n2 = mol.GetAtomWithIdx(match[2])
        if sulfonyl_on_neighbor(n1, core_atoms):
            return True, ("Found a urea core where one of the nitrogen atoms is substituted (directly or via an intervening NH) "
                          "by a sulfonyl group S(=O)(=O) with proper bonding")
        if sulfonyl_on_neighbor(n2, core_atoms):
            return True, ("Found a urea core where one of the nitrogen atoms is substituted (directly or via an intervening NH) "
                          "by a sulfonyl group S(=O)(=O) with proper bonding")
    
    return False, "Urea core found, but no nitrogen is substituted by a sulfonyl group (S(=O)(=O)-R or NH-S(=O)(=O)-R)"

# For testing purposes:
if __name__ == "__main__":
    test_smiles = [
        # chlorsulfuron: previously a false negative; now should be recognized.
        "C1(=NC(=NC(=N1)NC(NS(C=2C(=CC=CC2)Cl)(=O)=O)=O)OC)C",
        # 3-[(2-adamantylamino)-oxomethyl]-1-methylsulfonyl-1-pentylurea: true positive.
        "CCCCCN(C(=O)NC(=O)NC1C2CC3CC1CC(C2)C3)S(=O)(=O)C",
        # glipizide: true positive.
        "Cc1cnc(cn1)C(=O)NCCc1ccc(cc1)S(=O)(=O)NC(=O)NC1CCCCC1",
        # BM 567 (false positive before) – should not classify as N-sulfonylurea.
        "S(=O)(=O)(NC(=O)NCCCCC)C1=C(NC2CCCCC2)C=CC([N+]([O-])=O)=C1",
    ]
    for s in test_smiles:
        res, reason = is_N_sulfonylurea(s)
        print(f"SMILES: {s}\nResult: {res}\nReason: {reason}\n")