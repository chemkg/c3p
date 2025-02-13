"""
Classifies: CHEBI:46761 dipeptide
"""
"""
Classifies: dipeptide
Definition: A dipeptide is any molecule that contains exactly two amino‐acid residues
connected (directly or in a cyclic fashion) by (one or two) peptide linkage(s).
This implementation finds candidate amide bonds (where a carbonyl carbon is bonded to a nitrogen)
and then extracts the alpha–carbons (the alpha carbon being a carbon attached to the carbonyl in one residue
or attached to the nitrogen in the other). Finally, it verifies that the two candidate alpha carbons 
are connected by a 3–bond path: alpha – carbonyl – nitrogen – alpha.
If exactly one unique pair is found the function reports a positive dipeptide classification.
Otherwise it gives a negative classification.
Note: This heuristic may miss some valid dipeptides or falsely include non‐dipeptides.
"""

from rdkit import Chem
from rdkit.Chem import rdMolOps

def is_dipeptide(smiles: str):
    """
    Determines if a molecule is a dipeptide based on its SMILES string.
    A dipeptide is defined as any molecule that contains two amino‐acid residues connected by peptide linkage(s).

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if molecule is classified as a dipeptide, False otherwise.
        str: Reason for classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    candidate_pairs = set()  # each element will be a frozenset of two atom indices (the candidate alpha carbons)
    
    # Loop over all bonds and look for peptide bond candidates:
    for bond in mol.GetBonds():
        # We want a single bond connecting a carbon and a nitrogen
        if bond.GetBondType() != Chem.rdchem.BondType.SINGLE:
            continue
        
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        # Identify the carbon and nitrogen in the bond.
        if a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 7:
            carbon_atom = a1
            nitrogen_atom = a2
        elif a1.GetAtomicNum() == 7 and a2.GetAtomicNum() == 6:
            nitrogen_atom = a1
            carbon_atom = a2
        else:
            continue
        
        # Check that the carbon is "carbonyl-like": it must have at least one double bond to an oxygen.
        has_carbonyl = False
        for bond_c in carbon_atom.GetBonds():
            if bond_c.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                other = bond_c.GetOtherAtom(carbon_atom)
                if other.GetAtomicNum() == 8:
                    has_carbonyl = True
                    break
        if not has_carbonyl:
            continue
        
        # Now try to get the candidate alpha carbon from the carbonyl side.
        alpha_from_C = None
        for nb in carbon_atom.GetNeighbors():
            # Exclude the nitrogen we are using as part of the peptide link and any oxygen (the carbonyl O)
            if nb.GetIdx() == nitrogen_atom.GetIdx():
                continue
            if nb.GetAtomicNum() == 6:
                alpha_from_C = nb.GetIdx()
                break
        if alpha_from_C is None:
            continue
        
        # For the nitrogen side, look for an attached carbon (other than the carbonyl carbon)
        alpha_from_N = None
        for nb in nitrogen_atom.GetNeighbors():
            if nb.GetIdx() == carbon_atom.GetIdx():
                continue
            if nb.GetAtomicNum() == 6:
                alpha_from_N = nb.GetIdx()
                break
        if alpha_from_N is None:
            continue
        
        # Now we have a candidate peptide bond linking residue with alpha carbon alpha_from_C
        # to residue with alpha carbon alpha_from_N.
        # In a proper peptide linkage these two alpha carbons should be separated by exactly three bonds:
        # (alpha_from_C) -- (carbonyl C) -- (amide N) -- (alpha_from_N)
        path = rdMolOps.GetShortestPath(mol, alpha_from_C, alpha_from_N)
        if len(path) - 1 != 3:
            # If the bond path length is not exactly 3, skip this candidate.
            continue
        
        # Add the unique pair (order does not matter).
        candidate_pairs.add(frozenset({alpha_from_C, alpha_from_N}))
    
    # For dipeptides, we expect a single unique peptide linkage connecting two amino acid residues.
    if not candidate_pairs:
        return False, "No peptide linkage connecting two amino–acid residues found"
    elif len(candidate_pairs) == 1:
        return True, "Found two amino–acid residues connected by a peptide linkage"
    else:
        return False, f"Found multiple distinct peptide linkages: {len(candidate_pairs)}; not a dipeptide"

# Example usage (to test, uncomment):
# test_smiles = "SC[C@H](NC(=O)[C@@H](N)CCCCN)C(O)=O"  # Lys-Cys example
# result, reason = is_dipeptide(test_smiles)
# print(result, reason)

# Note: This heuristic may still fail for modified dipeptides where extra alpha–carbon candidates are present
# or in cases of extended backbone modifications. Further refinement (or even a machine–learning approach)
# could improve the classification.