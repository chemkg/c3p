"""
Classifies: CHEBI:50760 N-hydroxy-alpha-amino-acid
"""
"""
Classifies: N-hydroxy-alpha-amino-acid 
Definition: Any amino acid in which at least one hydrogen attached to an amino group 
           (either the backbone amino or a side-chain amine of an amino acid) is replaced by a hydroxy group.
           
Approach:
  1. Confirm that the molecule has an amino acid backbone. We do this by matching a SMARTS pattern
     for an α–carbon connected to both an amino (N) group and a carboxyl group. Our pattern tolerates
     carboxylate (deprotonated) or acid forms.
  2. Verify that the molecule contains at least one N–OH group. To do this, after adding hydrogens,
     we iterate over all nitrogen atoms and check whether any one has a neighboring oxygen via a single bond
     that in turn carries at least one hydrogen.
     
If both conditions are met, we classify the input as N–hydroxy–α–amino–acid.
"""

from rdkit import Chem

def is_N_hydroxy_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule (by SMILES) is an N-hydroxy-α-amino-acid.
    
    We require:
      1. The molecule contains an amino acid backbone: an α–carbon that is connected to both an amino group and
         a carboxyl group. Carboxyl groups are recognized as C(=O)[O] or C(=O)[O-].
      2. At least one nitrogen (either the backbone amine or a side-chain amine) carries a hydroxy substituent,
         i.e. it is covalently attached (by a single bond) to an oxygen that itself bears at least one hydrogen.
    
    Args:
        smiles (str): SMILES string.
    
    Returns:
        bool: True if molecule qualifies as an N-hydroxy-α-amino-acid, False otherwise.
        str: Explanation.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # Work with explicit hydrogens to reliably check for -OH groups.
    mol = Chem.AddHs(mol)
    
    # --- Step 1: Check for an amino acid backbone ---
    # The SMARTS pattern looks for a carbon (alpha carbon) that is bonded to:
    #   - an amino nitrogen (any nitrogen, [NX3])
    #   - a carboxyl carbon which is defined as a carbon double-bonded to an oxygen and attached to an oxygen (which might be protonated or deprotonated)
    # This is intentionally simple and may not capture every edge case but works for standard amino acids.
    # The carboxyl group is represented as: C(=O)[O;H1,-] where [O;H1,-] matches oxygen with a hydrogen or a negative charge.
    aa_backbone_smarts = "[C]([NX3])(C(=O)[O;H1,-])"
    aa_backbone = Chem.MolFromSmarts(aa_backbone_smarts)
    if not mol.HasSubstructMatch(aa_backbone):
        return False, "No amino acid backbone detected (missing carbon linked to both an amine and a carboxyl group)."
    
    # --- Step 2: Look for an N–OH substitution ---
    def has_N_hydroxy(atom):
        # For a nitrogen atom, we look for at least one neighboring oxygen connected via a single bond that itself has at least one hydrogen.
        if atom.GetAtomicNum() != 7:
            return False
        for nbr in atom.GetNeighbors():
            # Consider only oxygen neighbors.
            if nbr.GetAtomicNum() != 8:
                continue
            bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
            if bond is None or bond.GetBondType() != Chem.BondType.SINGLE:
                continue
            # Check that this oxygen bears at least one hydrogen.
            # (After AddHs, explicit H neighbors are present.)
            # We also ensure that the oxygen is not double-bonded to a carbon (i.e. not a carbonyl oxygen) by checking its degree.
            h_count = sum(1 for Hb in nbr.GetNeighbors() if Hb.GetAtomicNum() == 1)
            if h_count > 0:
                return True
        return False

    # Now, iterate over all nitrogen atoms in the molecule.
    found_Nhydroxy = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7 and has_N_hydroxy(atom):
            found_Nhydroxy = True
            break

    if not found_Nhydroxy:
        return False, "No N–OH substituent detected on any nitrogen atom."
    
    return True, "Amino acid backbone detected and at least one nitrogen carries a hydroxy substituent."

    
# Example usage with the referenced SMILES examples:
if __name__ == "__main__":
    test_smiles = [
        "C(N)(=NO)NCCC[C@H](N)C(=O)O",        # N(5)-[amino(hydroxyimino)methyl]-L-ornithine
        "CSCCCCC(N(O)O)C(O)=O",               # N,N-dihydroxydihomomethionine
        "O=C(O)[C@@H](NO)CCCCCSC",            # N-hydroxy-L-trihomomethionine
        "N[C@@H](CCCCNO)C(O)=O",              # N(6)-hydroxy-L-lysine
        "CC[C@H](C)[C@H](N(O)O)C(O)=O",        # N,N-dihydroxy-L-isoleucine
        "O=C(O)[C@@H](NO)CCCCSC",             # N-hydroxy-L-dihomomethionine
        "ON(O)[C@@H](Cc1ccc(O)cc1)C(O)=O",     # N,N-dihydroxy-L-tyrosine
        "C(\\N)(=N/O)/NCCC[C@H](N)C(=O)O",     # N(5)-[(Z)-amino(hydroxyimino)methyl]-L-ornithine
        "NC(CCCN\\C(N)=N\\O)C(O)=O",           # N(5)-[(E)-amino(hydroxyimino)methyl]ornithine
        "O=N(C(O)=O)[C@@H](NO)CCCCCSC",        # Variant adjustment if needed
        "ONCC(O)=O",                         # N-hydroxyglycine
        "CC(C)[C@H](N(O)O)C(O)=O",             # N,N-dihydroxy-L-valine
        "CSCCCCCCCC(NO)C(O)=O",               # N-hydroxypentahomomethionine
        "O=C(O)[C@@H](NO)CCCCCCSC",            # N-hydroxy-L-tetrahomomethionine
        "ON(O)[C@@H](Cc1ccccc1)C(O)=O",        # N,N-dihydroxy-L-phenylalanine
        "CSCCCCCCCCC(NO)C(O)=O",              # N-hydroxyhexahomomethionine
        "CSCCCCCCCCC(N(O)O)C(O)=O",            # N,N-dihydroxyhexahomomethionine
        "CSCCCCCCCC(N(O)O)C(O)=O",             # N,N-dihydroxypentahomomethionine
        "NC(CCCNC(=N)NO)C(O)=O",               # N(5)-[(hydroxyamino)(imino)methyl]ornithine
        "N[C@@H](CCCNC(=N)NO)C(O)=O",          # N(5)-[(hydroxyamino)(imino)methyl]-L-ornithine
        "ON(O)[C@@H](Cc1c[nH]c2ccccc12)C(O)=O", # N,N-dihydroxy-L-tryptophan
        "CC(C)[C@H](NO)C(O)=O",               # N-hydroxy-L-valine
        "CSCCCCCC(N(O)O)C(O)=O",              # N,N-dihydroxytrihomomethionine
        "N1([C@@H](CCCC1)C(=O)O)O",            # N-hydroxy-L-pipecolic acid
    ]
    
    for smi in test_smiles:
        result, reason = is_N_hydroxy_alpha_amino_acid(smi)
        print("SMILES:", smi)
        print("Result:", result)
        print("Reason:", reason)
        print("--------------------------------------------------")