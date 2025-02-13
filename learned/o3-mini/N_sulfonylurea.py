"""
Classifies: CHEBI:76983 N-sulfonylurea
"""
"""
Classifies: N-sulfonylurea
A urea in which one of the hydrogens attached to a nitrogen of the urea group is replaced by a sulfonyl group.
This key moiety is found in various herbicides and antidiabetic drugs.
"""

from rdkit import Chem

def is_N_sulfonylurea(smiles: str):
    """
    Determines if a molecule is a N-sulfonylurea based on its SMILES string.
    A N-sulfonylurea is defined as a urea (N-C(=O)-N) where one of the urea nitrogens is substituted (directly or
    via a single intervening amino group) by a sulfonyl group (-S(=O)(=O)-R).
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is classified as N-sulfonylurea, False otherwise.
        str: Explanation (reason) for the classification.
    """
    # Parse the input SMILES.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Use a simplified urea pattern that finds a urea fragment.
    # This pattern “N-C(=O)-N” should match even if extra substituents are present.
    urea_smarts = "N-C(=O)-N"
    urea_core = Chem.MolFromSmarts(urea_smarts)
    if urea_core is None:
        return False, "Failed to create urea SMARTS query"
    
    urea_matches = mol.GetSubstructMatches(urea_core)
    if not urea_matches:
        return False, "No urea (N-C(=O)-N) core found in the molecule"
    
    def is_sulfonyl(sulfur_atom):
        """
        Checks whether a sulfur atom looks like a sulfonyl center.
        We require two oxygens connected by double bonds and that the sulfur has a total of 4 connections.
        """
        if sulfur_atom.GetAtomicNum() != 16:
            return False
        dbl_oxygen_count = 0
        for bond in sulfur_atom.GetBonds():
            nbr = bond.GetOtherAtom(sulfur_atom)
            if nbr.GetAtomicNum() == 8 and bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                dbl_oxygen_count += 1
        # A proper sulfonyl group characteristically has two S=O double bonds and total degree 4.
        if dbl_oxygen_count == 2 and sulfur_atom.GetDegree() == 4:
            return True
        return False

    def sulfonyl_on_neighbor(n_atom, core_atoms, visited=None):
        """
        Checks whether the given nitrogen (n_atom) (which is part of the urea core) has a substituent that
        eventually leads (within one intervening amino) to a sulfonyl S.
        
        Args:
            n_atom: the urea nitrogen atom
            core_atoms: indices of atoms in the urea core (so that we ignore intramolecular urea bonds)
            visited: set of visited atom indices to avoid loops
        
        Returns:
            True if a valid sulfonyl connection is found, False otherwise.
        """
        if visited is None:
            visited = set()
        visited.add(n_atom.GetIdx())
        for nbr in n_atom.GetNeighbors():
            # Skip if neighbor is part of the urea core.
            if nbr.GetIdx() in core_atoms:
                continue
            # CASE 1: Direct attachment to a sulfur.
            if nbr.GetAtomicNum() == 16:
                bond = mol.GetBondBetweenAtoms(n_atom.GetIdx(), nbr.GetIdx())
                if bond is not None and bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                    if is_sulfonyl(nbr):
                        return True
            # CASE 2: One intervening amino (nitrogen) group.
            if nbr.GetAtomicNum() == 7 and nbr.GetIdx() not in visited:
                bond = mol.GetBondBetweenAtoms(n_atom.GetIdx(), nbr.GetIdx())
                if bond is None or bond.GetBondType() != Chem.rdchem.BondType.SINGLE:
                    continue
                # To be considered a simple intervening amino, limit its heavy-atom connectivity.
                # We require that it has exactly 2 heavy-atom neighbors (the urea N and the sulfonyl S).
                heavy_neighbors = [a for a in nbr.GetNeighbors() if a.GetAtomicNum() > 1]
                if len(heavy_neighbors) > 2:
                    # Overly substituted intervening amino; skip it.
                    continue
                # Now look at neighbors (other than n_atom and core atoms) from the intervening nitrogen.
                for nbr2 in nbr.GetNeighbors():
                    if nbr2.GetIdx() == n_atom.GetIdx() or nbr2.GetIdx() in core_atoms:
                        continue
                    if nbr2.GetAtomicNum() == 16:
                        bond2 = mol.GetBondBetweenAtoms(nbr.GetIdx(), nbr2.GetIdx())
                        if bond2 is not None and bond2.GetBondType() == Chem.rdchem.BondType.SINGLE:
                            if is_sulfonyl(nbr2):
                                return True
        return False

    # Look at every match for the urea fragment.
    # The simple SMARTS "N-C(=O)-N" returns a 3-atom match: first N, the carbonyl C, and second N.
    for match in urea_matches:
        core_atoms = set(match)
        # Grab the two terminal nitrogens.
        n1 = mol.GetAtomWithIdx(match[0])
        n2 = mol.GetAtomWithIdx(match[2])
        # Check if either terminal nitrogen carries (directly or via a simple NH linker) a sulfonyl group.
        if sulfonyl_on_neighbor(n1, core_atoms):
            return True, ("Found a urea core with a terminal nitrogen that is substituted (directly or via a simple NH linker) "
                          "by a sulfonyl group S(=O)(=O)-R")
        if sulfonyl_on_neighbor(n2, core_atoms):
            return True, ("Found a urea core with a terminal nitrogen that is substituted (directly or via a simple NH linker) "
                          "by a sulfonyl group S(=O)(=O)-R")
    
    return False, ("Urea core found, but neither urea nitrogen is substituted (directly or via a simple NH linker) "
                   "by a sulfonyl group S(=O)(=O)-R")

# For testing purposes:
if __name__ == "__main__":
    test_smiles = [
        # chlorsulfuron (previously false negative, should now be recognized)
        "C1(=NC(=NC(=N1)NC(NS(C=2C(=CC=CC2)Cl)(=O)=O)=O)OC)C",
        # True positive: 3-[(2-adamantylamino)-oxomethyl]-1-methylsulfonyl-1-pentylurea
        "CCCCCN(C(=O)NC(=O)NC1C2CC3CC1CC(C2)C3)S(=O)(=O)C",
        # True positive: glipizide
        "Cc1cnc(cn1)C(=O)NCCc1ccc(cc1)S(=O)(=O)NC(=O)NC1CCCCC1",
        # False positive test (BM 567) – should not be classified as N-sulfonylurea.
        "S(=O)(=O)(NC(=O)NCCCCC)C1=C(NC2CCCCC2)C=CC([N+]([O-])=O)=C1",
    ]
    for s in test_smiles:
        res, reason = is_N_sulfonylurea(s)
        print(f"SMILES: {s}\nResult: {res}\nReason: {reason}\n")