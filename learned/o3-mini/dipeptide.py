"""
Classifies: CHEBI:46761 dipeptide
"""
"""
Classifies: dipeptide
Definition: A dipeptide is any molecule that contains two amino‐acid residues connected
(directly or in a cyclic fashion) by peptide linkage(s).

The heuristic implemented here does the following:
 - Parse the SMILES string.
 - Loop through all bonds in the molecule and look for a single C–N bond.
 - For such bonds, ensure the carbon is “carbonyl-like” (has a double bond to oxygen).
 - Identify a candidate α–carbon attached to the carbon (excluding the carbonyl O and the N)
   and a candidate α–carbon attached to the nitrogen (excluding the carbonyl C).
 - Use the RDKit shortest‐path function (via rdmolops) to verify that the two candidate α–carbons 
   are separated by three bonds (as in the ideal linear peptide bond arrangement).
 - If exactly one unique candidate pair is found, classify as a dipeptide.
Note: This heuristic is simple and may miss some valid structures (especially in cyclic cases)
or wrongly classify modified dipeptides.
"""

from rdkit import Chem
from rdkit.Chem import rdmolops

def is_dipeptide(smiles: str):
    """
    Determines if a molecule is a dipeptide based on its SMILES string.
    A dipeptide is defined as any molecule that contains exactly two amino‐acid residues connected by peptide linkage(s).

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if molecule is classified as a dipeptide, False otherwise.
        str: Reason for classification.
    """
    # Parse SMILES into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    candidate_pairs = set()  # each element is a frozenset of two candidate alpha-carbon atom indices
    
    # Loop over all bonds to find candidate peptide bonds
    for bond in mol.GetBonds():
        # We work with a single bond connecting a carbon and a nitrogen
        if bond.GetBondType() != Chem.rdchem.BondType.SINGLE:
            continue

        atom1 = bond.GetBeginAtom()
        atom2 = bond.GetEndAtom()
        
        # Determine which atom is carbon and which is nitrogen
        if atom1.GetAtomicNum() == 6 and atom2.GetAtomicNum() == 7:
            carbon_atom = atom1
            nitrogen_atom = atom2
        elif atom1.GetAtomicNum() == 7 and atom2.GetAtomicNum() == 6:
            nitrogen_atom = atom1
            carbon_atom = atom2
        else:
            continue

        # Check that the carbon is carbonyl-like (has a double bond to an oxygen)
        has_carbonyl = False
        for cbond in carbon_atom.GetBonds():
            if cbond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                other_at = cbond.GetOtherAtom(carbon_atom)
                if other_at.GetAtomicNum() == 8:
                    has_carbonyl = True
                    break
        if not has_carbonyl:
            continue

        # Look for candidate alpha-carbon on the carbonyl side:
        # It should be a carbon neighbor of the carbon excluding the peptide nitrogen and the oxygen (from the carbonyl).
        alpha_from_C = None
        for neighbor in carbon_atom.GetNeighbors():
            if neighbor.GetIdx() == nitrogen_atom.GetIdx():
                continue
            if neighbor.GetAtomicNum() == 6:
                alpha_from_C = neighbor.GetIdx()
                break
        if alpha_from_C is None:
            continue

        # For the nitrogen side, look for an attached carbon other than the carbonyl carbon.
        alpha_from_N = None
        for neighbor in nitrogen_atom.GetNeighbors():
            if neighbor.GetIdx() == carbon_atom.GetIdx():
                continue
            if neighbor.GetAtomicNum() == 6:
                alpha_from_N = neighbor.GetIdx()
                break
        if alpha_from_N is None:
            continue

        # Verify that the candidate alpha-carbons are connected by a 3-bond path:
        # Ideal linear dipeptide gives: α-carbon – carbonyl C – amide N – α-carbon.
        try:
            path = rdmolops.GetShortestPath(mol, alpha_from_C, alpha_from_N)
        except Exception as e:
            return False, f"Error computing shortest path: {e}"
        # The number of bonds is (number of atoms in path - 1)
        if len(path) - 1 != 3:
            continue

        # Record the unique candidate pair (order does not matter)
        candidate_pairs.add(frozenset({alpha_from_C, alpha_from_N}))

    # Evaluate results based on the number of unique peptide linkage candidate pairs.
    if not candidate_pairs:
        return False, "No peptide linkage connecting two amino–acid residues found"
    elif len(candidate_pairs) == 1:
        return True, "Found two amino–acid residues connected by a peptide linkage"
    else:
        return False, f"Found multiple distinct peptide linkages ({len(candidate_pairs)}); not a unique dipeptide linkage"

# Example usage for testing:
if __name__ == "__main__":
    # Try one dipeptide example, e.g., Lys-Cys:
    test_smiles = "SC[C@H](NC(=O)[C@@H](N)CCCCN)C(O)=O"
    result, reason = is_dipeptide(test_smiles)
    print(result, reason)