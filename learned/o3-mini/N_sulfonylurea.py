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
    by a sulfonyl group (–S(=O)2R).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is N-sulfonylurea, False otherwise
        str: Reason for the classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS for a simple urea core: N-C(=O)-N.
    # Note: This pattern may match other urea derivatives as well.
    urea_smarts = "NC(=O)N"
    urea_core = Chem.MolFromSmarts(urea_smarts)
    if urea_core is None:
        return False, "Invalid urea SMARTS pattern"
    
    # Find urea substructure matches.
    urea_matches = mol.GetSubstructMatches(urea_core)
    if not urea_matches:
        return False, "No urea (N-C(=O)-N) core found in the molecule"
    
    # For each found urea core, check if one of the N atoms is substituted with a sulfonyl group.
    # We'll inspect neighbors of the urea N atoms (the first and third atoms in the match)
    for match in urea_matches:
        # In the SMARTS "NC(=O)N", match indices:
        # match[0] -> first N
        # match[1] -> carbonyl carbon
        # match[2] -> second N
        for n_idx in (match[0], match[2]):
            n_atom = mol.GetAtomWithIdx(n_idx)
            # Check each neighbor of this nitrogen
            for neighbor in n_atom.GetNeighbors():
                # Skip the carbonyl carbon already in the urea core
                if neighbor.GetAtomicNum() == 6:
                    continue
                # Look for a sulfur neighbor (atomic number 16) attached via a single bond
                if neighbor.GetAtomicNum() == 16:
                    # Check if the bond between the nitrogen and sulfur is a single bond.
                    bond = mol.GetBondBetweenAtoms(n_atom.GetIdx(), neighbor.GetIdx())
                    if bond is None or bond.GetBondType() != Chem.rdchem.BondType.SINGLE:
                        continue
                    # Now check that this sulfur has at least two double-bonded oxygens.
                    sulfonyl_O_count = 0
                    for sulfon_neighbor in neighbor.GetNeighbors():
                        # if the neighbor is an oxygen and the bond is double
                        if sulfon_neighbor.GetAtomicNum() == 8:
                            s_bond = mol.GetBondBetweenAtoms(neighbor.GetIdx(), sulfon_neighbor.GetIdx())
                            if s_bond is not None and s_bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                                sulfonyl_O_count += 1
                    if sulfonyl_O_count >= 2:
                        return True, ("Found a urea core with a nitrogen substituted by a sulfonyl group "
                                      "(has sulfur with at least two double-bonded oxygens)")
    
    # If no suitable substitution found, then this is not an N-sulfonylurea.
    return False, "Urea core found, but no nitrogen has a sulfonyl substitution (S(=O)(=O)-R) attached"

# Below are some examples for testing the function.
if __name__ == "__main__":
    test_smiles = [
        # Chlorsulfuron
        "C1(=NC(=NC(=N1)NC(NS(C=2C(=CC=CC2)Cl)(=O)=O)=O)OC)C",
        # 3-[(2-adamantylamino)-oxomethyl]-1-methylsulfonyl-1-pentylurea
        "CCCCCN(C(=O)NC(=O)NC1C2CC3CC1CC(C2)C3)S(=O)(=O)C",
        # glipizide
        "Cc1cnc(cn1)C(=O)NCCc1ccc(cc1)S(=O)(=O)NC(=O)NC1CCCCC1",
    ]
    for s in test_smiles:
        result, reason = is_N_sulfonylurea(s)
        print(f"SMILES: {s}\nResult: {result}\nReason: {reason}\n")