"""
Classifies: CHEBI:46761 dipeptide
"""
"""
Classifies: dipeptide
Definition: A dipeptide is any molecule that contains two amino‐acid residues 
connected (directly or in a cyclic fashion) by peptide linkage(s).
This implementation uses a heuristic based on identifying candidate peptide bonds:
  - It scans for C–N single bonds where the carbon shows a carbonyl (C=O) bond.
  - For such bonds, it identifies a candidate α–carbon on the carbon side (excluding the 
    attached nitrogen and oxygen) and a candidate α–carbon on the nitrogen side (excluding the 
    carbonyl carbon).
  - It then verifies that these two candidate α–carbons are connected by a 3-bond path 
    (i.e. α–carbon – carbonyl C – amide N – α–carbon, as expected for a typical peptide bond).
  - Rather than forcing exactly one candidate peptide bond, it collects all candidate α–carbon 
    indices. A valid dipeptide should include exactly two distinct α–carbons (even if the peptide 
    linkage is redundant, as in cyclic dipeptides).
The function returns True if exactly two distinct α–carbons are found and if at least one connecting 
peptide bond is detected (with a 3–bond path between the α–carbons). Otherwise it returns False.
"""

from rdkit import Chem
from rdkit.Chem import rdmolops

def is_dipeptide(smiles: str):
    """
    Determines if a molecule is a dipeptide based on its SMILES string.
    A dipeptide is defined as any molecule that contains two amino–acid residues connected by peptide linkage(s).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a dipeptide, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # We will collect candidate peptide bond pairs.
    # candidate_pairs will be a list of tuples (alpha_from_C, alpha_from_N) of candidate α–carbon atom indices.
    candidate_pairs = []
    # Also gather a set of all candidate α–carbon indices.
    alpha_set = set()

    # Loop over all bonds in the molecule looking for candidate peptide bonds.
    for bond in mol.GetBonds():
        # Consider only single bonds
        if bond.GetBondType() != Chem.rdchem.BondType.SINGLE:
            continue

        atom1 = bond.GetBeginAtom()
        atom2 = bond.GetEndAtom()

        # We want a bond between a carbon and a nitrogen.
        if atom1.GetAtomicNum() == 6 and atom2.GetAtomicNum() == 7:
            carbon_atom = atom1
            nitrogen_atom = atom2
        elif atom1.GetAtomicNum() == 7 and atom2.GetAtomicNum() == 6:
            nitrogen_atom = atom1
            carbon_atom = atom2
        else:
            continue

        # Check that the carbon atom is "carbonyl-like": has at least one double bond to oxygen.
        has_carbonyl = False
        for cbond in carbon_atom.GetBonds():
            if cbond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                other_atom = cbond.GetOtherAtom(carbon_atom)
                if other_atom.GetAtomicNum() == 8:
                    has_carbonyl = True
                    break
        if not has_carbonyl:
            continue

        # Identify candidate α–carbon on the carbon side.
        # Choose a neighboring carbon of the carbonyl carbon that is not the peptide nitrogen
        # and not the oxygen from the carbonyl.
        alpha_from_C = None
        for nbr in carbon_atom.GetNeighbors():
            # Skip if this neighbor is the peptide nitrogen.
            if nbr.GetIdx() == nitrogen_atom.GetIdx():
                continue
            # Also skip if the neighbor is oxygen (likely the carbonyl O).
            if nbr.GetAtomicNum() == 8:
                continue
            if nbr.GetAtomicNum() == 6:
                alpha_from_C = nbr.GetIdx()
                break
        if alpha_from_C is None:
            continue

        # Identify candidate α–carbon on the nitrogen side.
        # Look for a neighboring carbon (other than the carbonyl carbon) on the nitrogen side.
        alpha_from_N = None
        for nbr in nitrogen_atom.GetNeighbors():
            if nbr.GetIdx() == carbon_atom.GetIdx():
                continue
            if nbr.GetAtomicNum() == 6:
                alpha_from_N = nbr.GetIdx()
                break
        if alpha_from_N is None:
            continue

        # For a candidate peptide bond, check that the two candidate α–carbons are connected by a 3–bond path.
        try:
            path = rdmolops.GetShortestPath(mol, alpha_from_C, alpha_from_N)
        except Exception as e:
            return False, f"Error computing shortest path: {e}"
        # A typical peptide bond (linear case) gives a path: α–C -> carbonyl C -> amide N -> α–N.
        if len(path) - 1 != 3:
            continue

        # Record this candidate peptide bond (the order of α–carbons does not matter).
        candidate_pairs.append( (alpha_from_C, alpha_from_N) )
        alpha_set.update([alpha_from_C, alpha_from_N])
    
    # Evaluate: For a dipeptide we expect that the union of all candidate α–carbons equals two.
    if len(candidate_pairs) == 0:
        return False, "No peptide linkage connecting two amino–acid residues found"
    
    if len(alpha_set) == 2:
        # Optionally, count how many candidate bonds join these same two residues.
        num_candidate_bonds = sum(1 for pair in candidate_pairs if set(pair) == alpha_set)
        # For a linear dipeptide num_candidate_bonds should be 1 and for a cyclic dipeptide it may be 2.
        if num_candidate_bonds in [1, 2]:
            return True, "Found two amino–acid residues connected by peptide linkage(s)"
        else:
            return False, f"Found {num_candidate_bonds} connecting bonds between two residues; unusual pattern."
    else:
        return False, f"Found multiple distinct amino–acid residues ({len(alpha_set)}); not a unique dipeptide linkage"

# Example usage for testing:
if __name__ == "__main__":
    # Test with a valid dipeptide example, e.g., Lys-Cys.
    test_smiles = "SC[C@H](NC(=O)[C@@H](N)CCCCN)C(O)=O"
    result, reason = is_dipeptide(test_smiles)
    print(result, reason)
    
    # You can test additional examples as needed.
    
"""
Note:
This heuristic does not catch every possible dipeptide (or modified / cyclic variant)
and may also mis‐classify a few borderline cases. Adjustments (or a more thorough fragmentation
approach) might be needed if further refinement is desired.
"""