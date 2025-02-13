"""
Classifies: CHEBI:16385 organic sulfide
"""
"""
Classifies: Organic sulfide (thioether) compounds.
Definition: Compounds having the structure R–S–R (with R ≠ H).
Here we require that an S atom (atomic number 16) in the molecule is attached
to exactly two heavy atoms, both of which are carbon (atomic number 6). In addition,
we “veto” molecules that show a clear peptide backbone (by detecting at least two
amide bonds), which are common false‐positives (e.g. methionine residues in peptides).
Note: This is a heuristic implementation and may not capture all edge cases.
"""

from rdkit import Chem

def is_organic_sulfide(smiles: str):
    """
    Determines if a molecule is an organic sulfide (thioether) based on its SMILES string.
    Here the criterion is that the molecule contains at least one sulfur atom that is bonded
    to exactly two non-hydrogen atoms and, for improved discrimination, both neighbors must be carbons.
    Also, if the molecule contains multiple peptide-bond motifs, it is likely a peptide and we do not
    wish to classify that as “organic sulfide” (even if a thioether appears in a side‐chain).
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as an organic sulfide, False otherwise.
        str: A reason string describing the decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Crude peptide detection: count amide bonds (C(=O)N).
    amide_smarts = Chem.MolFromSmarts("C(=O)N")
    peptide_matches = mol.GetSubstructMatches(amide_smarts) if amide_smarts is not None else []
    peptide_flag = len(peptide_matches) > 1  # if >1 amide bond exists, we treat it as peptide-like
    
    # Loop over all atoms; look for a thioether sulfur:
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 16:
            # Check number of neighbors (here we work on the implicit (not added) molecule).
            neighbors = atom.GetNeighbors()
            # Expect exactly 2 bonds for an isolated thioether sulfur.
            if len(neighbors) != 2:
                continue
            # Check none of the neighbors is a hydrogen (atomic num 1)
            if any(nbr.GetAtomicNum() == 1 for nbr in neighbors):
                continue
            # To be sure the S is in a typical organic sulfide motif, require that both neighbors are carbon.
            if not all(nbr.GetAtomicNum() == 6 for nbr in neighbors):
                continue
            # Also, if the molecule strongly appears to be a peptide (e.g. several amide bonds),
            # we decide not to classify it as an organic sulfide.
            if peptide_flag:
                # We could refine by checking if the S is within a “side-chain” of a peptide,
                # but here we take a simple stance.
                continue
            # If all conditions are met, we classify as an organic sulfide.
            return True, "Contains an R–S–R motif (S bonded to two carbons, with no S–H bonds) and no strong peptide signals"
    
    # If we get here no sulfur in the proper environment was found.
    return False, "No R–S–R (thioether) motif found"

# Example usage (uncomment to run some tests):
# test_smiles = [
#     "CSC",                      # dimethyl sulfide, should be True
#     "SC",                       # methanethiol, should be False (has S–H bond)
#     "c1ccsc1",                  # thiophene, should be True
#     "S(CCCCC)C",                # 1-(Methylthio)pentane, should be True
#     "N[C@@H](CSC(C#N)c1c[nH]c2ccccc12)C(O)=O",  # L-Cys(IAN), should be True
#     "C1=CC=C(C=C1)CCCNC2=[NH+]C=NC3=C2NC=N3"      # LSM-21130*, no S so will be False
# ]
#
# for s in test_smiles:
#     result, reason = is_organic_sulfide(s)
#     print(f"SMILES: {s}\n  -> {result}: {reason}\n")