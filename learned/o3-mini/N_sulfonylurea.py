"""
Classifies: CHEBI:76983 N-sulfonylurea
"""
"""
Classifies: N-sulfonylurea
A urea in which one of the hydrogens attached to a nitrogen of the urea group is replaced by a sulfonyl group.
This moiety is common in various herbicides and antidiabetic drugs.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

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
        str: Reason for classification
    """
    # Parse the input SMILES into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # First search for a urea core. Use a SMARTS that matches N-C(=O)-N.
    # This is a minimal pattern for a urea group.
    urea_smarts = "N-C(=O)-N"
    urea_core = Chem.MolFromSmarts(urea_smarts)
    if urea_core is None:
        return False, "Could not interpret urea SMARTS pattern"
    
    urea_matches = mol.GetSubstructMatches(urea_core)
    if not urea_matches:
        return False, "No urea (N-C(=O)-N) core found in the molecule"
    
    def check_sulfonyl(sulfur_atom):
        """
        Check if a sulfur atom is a sulfonyl moiety by counting oxygen neighbors
        that are double-bonded.
        """
        double_oxygens = 0
        # Loop through each neighbor of the sulfur.
        for neigh in sulfur_atom.GetNeighbors():
            if neigh.GetAtomicNum() == 8:
                bond = mol.GetBondBetweenAtoms(sulfur_atom.GetIdx(), neigh.GetIdx())
                if bond is not None and bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                    double_oxygens += 1
        return double_oxygens >= 2

    def sulfonyl_on_neighbor(n_atom, visited=set()):
        """
        Check if the nitrogen atom (n_atom) has a substituent which is (or leads to) a sulfonyl group.
        We ignore the bond that is part of the urea core.
        To allow for cases where an intervening amino (N) group exists, we check one level deeper.
        visited: to avoid loops.
        """
        # Mark current nitrogen as visited.
        visited = visited | {n_atom.GetIdx()}
        for nbr in n_atom.GetNeighbors():
            # Skip if neighbor is a carbonyl carbon (C with =O bond to our urea core)
            if nbr.GetAtomicNum() == 6:
                continue
            # Direct S check
            if nbr.GetAtomicNum() == 16:
                bond = mol.GetBondBetweenAtoms(n_atom.GetIdx(), nbr.GetIdx())
                if bond is not None and bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                    if check_sulfonyl(nbr):
                        return True
            # Allow one intervening nitrogen: if neighbor is nitrogen and we haven't visited it,
            # check its neighbors for a sulfur.
            if nbr.GetAtomicNum() == 7 and nbr.GetIdx() not in visited:
                for nbr2 in nbr.GetNeighbors():
                    if nbr2.GetIdx() == n_atom.GetIdx():
                        continue
                    if nbr2.GetAtomicNum() == 16:
                        bond2 = mol.GetBondBetweenAtoms(nbr.GetIdx(), nbr2.GetIdx())
                        if bond2 is not None and bond2.GetBondType() == Chem.rdchem.BondType.SINGLE:
                            if check_sulfonyl(nbr2):
                                return True
        return False

    # For each matched urea core, check the two terminal nitrogen atoms.
    # In the SMARTS "N-C(=O)-N", the match tuple has indices:
    # match[0] -> first N, match[1] -> C (carbonyl), match[2] -> second N.
    for match in urea_matches:
        # Check first nitrogen
        n1 = mol.GetAtomWithIdx(match[0])
        if sulfonyl_on_neighbor(n1):
            return True, ("Found a urea core with a nitrogen substituted (directly or via an intervening NH) "
                          "by a sulfonyl group (S with at least two double-bonded O atoms)")
        # Check second nitrogen
        n2 = mol.GetAtomWithIdx(match[2])
        if sulfonyl_on_neighbor(n2):
            return True, ("Found a urea core with a nitrogen substituted (directly or via an intervening NH) "
                          "by a sulfonyl group (S with at least two double-bonded O atoms)")
    
    return False, "Urea core found, but no nitrogen has a sulfonyl substitution (S(=O)(=O)-R or NH-S(=O)(=O)-R) attached"

# The code below is for testing purposes.
if __name__ == "__main__":
    test_smiles = [
        # chlorsulfuron -- should be recognized (previously false negative)
        "C1(=NC(=NC(=N1)NC(NS(C=2C(=CC=CC2)Cl)(=O)=O)=O)OC)C",
        # 3-[(2-adamantylamino)-oxomethyl]-1-methylsulfonyl-1-pentylurea
        "CCCCCN(C(=O)NC(=O)NC1C2CC3CC1CC(C2)C3)S(=O)(=O)C",
        # glipizide
        "Cc1cnc(cn1)C(=O)NCCc1ccc(cc1)S(=O)(=O)NC(=O)NC1CCCCC1",
    ]
    for s in test_smiles:
        res, reason = is_N_sulfonylurea(s)
        print(f"SMILES: {s}\nResult: {res}\nReason: {reason}\n")