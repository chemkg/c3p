"""
Classifies: CHEBI:16385 organic sulfide
"""
"""
Classifies: Organic sulfide (thioether) compounds.
Definition: Compounds having the structure R–S–R (with R ≠ H).
We require that at least one sulfur atom is present that is attached to exactly two heavy atoms,
both of which are carbon. In addition, to reduce mis‐classification of molecules 
that contain many amide (peptide) bonds we veto molecules with more than one C(=O)N fragment.
Furthermore, if either carbon neighbor is activated by an exocyclic carbonyl (C=O) – i.e. one
that is not “inside” a ring containing the S atom – then we will not count that –S– candidate.
Note: This is a heuristic implementation and may not capture all edge cases.
"""

from rdkit import Chem

def is_organic_sulfide(smiles: str):
    """
    Determines whether a molecule (given as a SMILES string) is an organic sulfide (thioether).
    
    The criteria used here are:
      1. The molecule must contain at least one sulfur atom (atomic number 16) that bridges exactly
         two atoms.
      2. Both atoms attached to that S must be carbons (atomic number 6) and none may be a hydrogen.
      3. To avoid mis‐classification of peptide or related structures, if the molecule contains more than one
         C(=O)N (amide) fragment then we assume the sulfide is part of a peptide.
      4. In addition, if one of the carbon neighbors is “activated” by an exocyclic carbonyl group – meaning
         the carbon has a double-bonded oxygen that is not in the same ring as the candidate S and that carbon –
         then that candidate is vetoed.
         
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as an organic sulfide, False otherwise.
        str: A reason string describing how that decision was reached.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 1. Count amide (C(=O)N) fragments.
    amide_smarts = Chem.MolFromSmarts("C(=O)N")
    peptide_matches = mol.GetSubstructMatches(amide_smarts) if amide_smarts else []
    if len(peptide_matches) > 1:
        # We veto if the overall molecule is peptide‐like.
        peptide_flag = True
    else:
        peptide_flag = False

    # 2. Get ring information (each ring is a tuple of atom indices)
    rings = list(mol.GetRingInfo().AtomRings())

    # 3. Look over all sulfur atoms in the molecule.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 16:
            continue
        # Check that the S atom has exactly two neighbors.
        neighbors = atom.GetNeighbors()
        if len(neighbors) != 2:
            continue
        # Reject candidate if any neighbor is hydrogen.
        if any(nbr.GetAtomicNum() == 1 for nbr in neighbors):
            continue
        # Both neighbors must be carbon.
        if not all(nbr.GetAtomicNum() == 6 for nbr in neighbors):
            continue

        veto_due_to_carbonyl = False
        # For each carbon neighbor, check for an externally attached carbonyl.
        for nbr in neighbors:
            for nbr2 in nbr.GetNeighbors():
                # Skip the bond going back to the candidate S
                if nbr2.GetIdx() == atom.GetIdx():
                    continue
                # If this bond is a double bond and the other atom is oxygen:
                bond = mol.GetBondBetweenAtoms(nbr.GetIdx(), nbr2.GetIdx())
                if bond is None:
                    continue
                if bond.GetBondType() == Chem.BondType.DOUBLE and nbr2.GetAtomicNum() == 8:
                    # Determine whether the oxygen is in the same ring as both the sulfur candidate and the carbon.
                    in_same_ring = False
                    for ring in rings:
                        if (atom.GetIdx() in ring) and (nbr.GetIdx() in ring) and (nbr2.GetIdx() in ring):
                            in_same_ring = True
                            break
                    # If the oxygen is not in the same ring, then veto this candidate.
                    if not in_same_ring:
                        veto_due_to_carbonyl = True
                        break
            if veto_due_to_carbonyl:
                break

        if veto_due_to_carbonyl:
            # This candidate S is adjacent to an exocyclic carbonyl on one of its carbons.
            continue

        if peptide_flag:
            # Even if the S passes, if the molecule is peptide-like, we skip.
            continue

        # We found at least one candidate S that is attached to two carbons, not to any hydrogen, and
        # neither carbon is externally activated by a nearby C=O. We therefore classify the molecule as an organic sulfide.
        return True, "Contains an R–S–R motif (thioether with S bonded to two carbons, no S–H) and lacks interfering carbonyl/peptide signals"

    # If we find no candidate S that passes our checks, we classify the molecule as not an organic sulfide.
    return False, "No valid R–S–R (thioether) motif found"

# Example usage (for testing, uncomment below):
# test_smiles_list = [
#     "COC(=O)CCSCCO[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O[C@@H]2O[C@H](CO)[C@H](O)[C@H](O)[C@H]2O)[C@H]1NC(C)=O",  # true positive
#     "S(CC(N(CC)CC)=O)C1=CC=C(N=N1)Cl",  # azintamide, should be true (was false negative before)
#     "S1CC(=O)C(O)=C1C",  # 4-Hydroxy-5-methyl-3(2H)-thiophenone, should be true
#     "S(C1=C(O)C=CC(=C1)CC(=O)OC)C2=C(O)C=CC(=C2)CC(=O)OC",  # Penicillithiophenol B, should be false
#     "CSC",  # dimethyl sulfide, simple true case.
# ]
# for s in test_smiles_list:
#     result, reason = is_organic_sulfide(s)
#     print(f"SMILES: {s}\n  -> {result}: {reason}\n")