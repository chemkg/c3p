"""
Classifies: CHEBI:16385 organic sulfide
"""
"""
Classifies: Organic sulfide (thioether) compounds.
Definition: Compounds having the structure R–S–R (with R ≠ H).
We require that at least one sulfur atom is present that is attached to exactly two heavy atoms,
both of which are carbon. In addition, to reduce mis‐classification of molecules 
that contain many amide (peptide) bonds (e.g. methionine residues in peptides) we veto molecules
with more than one C(=O)N fragment.
Furthermore, if either carbon neighbor is strongly activated by a carbonyl (C=O) group, we do not count 
that –S– candidate as a “stand‐alone” organic sulfide.
Note: This is a heuristic implementation and may not capture all edge cases.
"""

from rdkit import Chem

def is_organic_sulfide(smiles: str):
    """
    Determines if a molecule is an organic sulfide (thioether) based on its SMILES string.

    The criteria are:
      1. The molecule must contain at least one sulfur atom bridging two carbons
         (i.e. an R–S–R motif, with no S–H bonds).
      2. To avoid mis‐classifying peptide molecules, if the molecule contains more than one amide bond (C(=O)N)
         we assume the sulfide is part of a peptide.
      3. We try to avoid cases where one of the carbon neighbors is adjacent to a strongly electron‐withdrawing 
         carbonyl group (which may indicate the sulfide function is simply part of a larger functional group).

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is classified as an organic sulfide, False otherwise.
        str: A reason string describing the decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 1. Look for amide bonds (C(=O)N) as a crude indicator of peptide connectivity.
    amide_smarts = Chem.MolFromSmarts("C(=O)N")
    peptide_matches = mol.GetSubstructMatches(amide_smarts) if amide_smarts is not None else []
    peptide_flag = len(peptide_matches) > 1  # if more than one C(=O)N found, we treat molecule as peptide-like

    # 2. Use a SMARTS that looks for a S atom (atomic num 16) between two carbons (atomic num 6)
    #    that are not just hydrogens.
    thioether_smarts = Chem.MolFromSmarts("[#6;!H0]-S-[#6;!H0]")
    thioether_matches = mol.GetSubstructMatches(thioether_smarts)
    
    # 3. Go over each match and refine:
    #    - Ensure that the S atom indeed has exactly 2 neighbors.
    #    - Veto if any adjacent carbon is directly double-bonded to oxygen (as in carbonyls)
    #    - If the molecule clearly appears peptide-like, then skip.
    for match in thioether_matches:
        # match is a tuple: (atom_idx1, atom_idx_s, atom_idx2)
        c1_idx, s_idx, c2_idx = match
        s_atom = mol.GetAtomWithIdx(s_idx)
        if s_atom.GetAtomicNum() != 16:
            continue
        # Check that the S atom has exactly two neighbors.
        if len(s_atom.GetNeighbors()) != 2:
            continue
        # (Double check that neither neighbor is hydrogen)
        if any(nbr.GetAtomicNum() == 1 for nbr in s_atom.GetNeighbors()):
            continue
        
        # For each of the two carbon neighbors, check if they are directly attached to an oxygen via a double bond.
        veto_due_to_carbonyl = False
        for nbr in s_atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6:
                for other in nbr.GetNeighbors():
                    if other.GetAtomicNum() == 8:
                        bond = mol.GetBondBetweenAtoms(nbr.GetIdx(), other.GetIdx())
                        if bond and bond.GetBondType() == Chem.BondType.DOUBLE:
                            veto_due_to_carbonyl = True
                            break
                if veto_due_to_carbonyl:
                    break
        if veto_due_to_carbonyl:
            continue

        # If the molecule appears peptide-like (many amide bonds), skip this candidate.
        if peptide_flag:
            continue

        # If we pass all checks, we found a valid thioether candidate.
        return True, "Contains an R–S–R motif (thioether with sulfur bonded to two carbons, no S–H) and lacks interfering carbonyl or peptide signals"

    # If no candidate passed the filters, we classify the molecule as not an organic sulfide.
    return False, "No valid R–S–R (thioether) motif found"

# Example usage (uncomment to test):
# test_smiles = [
#     "CSC",                      # dimethyl sulfide, expected True
#     "SC",                       # methanethiol, expected False (S–H bond present)
#     "c1ccsc1",                  # thiophene, expected True
#     "S(CCCCC)C",                # 1-(Methylthio)pentane, expected True
#     "N[C@@H](CSC(C#N)c1c[nH]c2ccccc12)C(O)=O",  # L-Cys(IAN), expected True
# ]
# for s in test_smiles:
#     result, reason = is_organic_sulfide(s)
#     print(f"SMILES: {s}\n  -> {result}: {reason}\n")