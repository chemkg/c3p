"""
Classifies: CHEBI:33709 amino acid
"""
"""
Classifies: Amino acids
Definition: A carboxylic acid containing one or more amino groups.

This improved approach:
  - Adds explicit hydrogens to help count amino hydrogens.
  - Checks for a carboxylic acid group (either protonated or deprotonated).
  - Iterates over nitrogen atoms to seek a “free” amino group (i.e. an N with at least one hydrogen 
    that is not exclusively bonded to carbonyl carbons).
  - Looks for an alpha amino acid backbone pattern using SMARTS.
  - Counts the number of amide bonds and alpha backbone atoms; if multiple such atoms (or many amide bonds)
    are found in a heavy molecule, the molecule is assumed to be a peptide.
    
Note: This heuristic may still mis‐classify edge cases.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_amino_acid(smiles: str):
    """
    Determines if a molecule can be classified as an amino acid based on its SMILES string.
    
    The molecule must have at least one carboxylic acid (–COOH or –COO–) group as well as one or more
    amino groups that are considered “free” (not part of an amide bond) or an isolated alpha backbone.
    If multiple alpha centers (or many amide linkages) are found and the weight is high, then the molecule 
    is likely to be a peptide.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as an amino acid, False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse the molecule from SMILES and add explicit hydrogens.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    mol = Chem.AddHs(mol)
    
    # --- Check for carboxylic acid group ---
    # Two SMARTS: one for protonated -COOH and one for deprotonated -COO-
    acid_smarts1 = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    acid_smarts2 = Chem.MolFromSmarts("[CX3](=O)[O-]")
    has_acid = mol.HasSubstructMatch(acid_smarts1) or mol.HasSubstructMatch(acid_smarts2)
    if not has_acid:
        return False, "No carboxylic acid group found"
    
    # --- Count amide bonds (which are typical for peptides) ---
    # Amide bond: a nitrogen single-bonded to a carbon that is double-bonded to an oxygen.
    amide_smarts = Chem.MolFromSmarts("[NX3][C](=O)")
    amide_matches = mol.GetSubstructMatches(amide_smarts)
    amide_count = len(amide_matches)
    
    # --- Look for an alpha amino acid backbone pattern ---
    # Pattern: a non-ring carbon bonded to an amino group and a carboxylic acid (protonated or deprotonated)
    alpha_patterns = [
        Chem.MolFromSmarts("[C;!R]([NX3])(C(=O)[O-])"),
        Chem.MolFromSmarts("[C;!R]([NX3])(C(=O)O)")
    ]
    alpha_matches = []
    alpha_found = False
    for pat in alpha_patterns:
        matches = mol.GetSubstructMatches(pat)
        if matches:
            alpha_matches.extend(matches)
    if len(alpha_matches) > 0:
        alpha_found = True
    
    # --- Detect free (non-amide) amino groups ---
    # We require that a nitrogen atom has at least one hydrogen and at least one bond 
    # that is not to a carbonyl carbon.
    free_amino_valid = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 7:
            continue
        # Count total hydrogens (explicit + implicit)
        if atom.GetTotalNumHs() < 1:
            continue
        # For each neighbor, check if at least one connection is not part of an amide bond.
        non_amide_connection = False
        for nbr in atom.GetNeighbors():
            # If neighbor is not carbon, we assume the connection is free.
            if nbr.GetAtomicNum() != 6:
                non_amide_connection = True
                break
            # For a carbon neighbor, check if that carbon is part of a carbonyl (C=O) bond.
            is_carbonyl = False
            for bond in nbr.GetBonds():
                # Skip the bond back to the nitrogen
                if bond.GetBeginAtomIdx() == atom.GetIdx() or bond.GetEndAtomIdx() == atom.GetIdx():
                    continue
                # Check if any bond from the carbon is a double bond to oxygen.
                other_atom = bond.GetOtherAtom(nbr)
                if other_atom.GetAtomicNum() == 8 and bond.GetBondTypeAsDouble() == 2.0:
                    is_carbonyl = True
                    break
            if not is_carbonyl:
                non_amide_connection = True
                break
        if non_amide_connection:
            free_amino_valid = True
            break

    # --- Peptide/Multiple alpha centers check ---
    # If more than one alpha backbone is found, or if there are multiple amide bonds in a high molecular weight molecule,
    # the structure is likely a peptide rather than an isolated amino acid.
    if (amide_count > 1 and not alpha_found):
        mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
        if mol_wt > 300:
            return False, "Molecule likely contains multiple peptide bonds without an isolated amino acid backbone"
    if len(alpha_matches) >= 2:
        return False, "Multiple alpha backbone atoms suggest a peptide rather than a single amino acid"
    
    # --- Combine results ---
    has_valid_amino = alpha_found or free_amino_valid
    if not has_valid_amino:
        return False, "No free (non-amide) amino group found"
    
    return True, "Molecule contains a carboxylic acid and a valid (free or alpha) amino group"

# Example usage (for testing):
if __name__ == "__main__":
    test_cases = [
        ("OC(=O)C(N)CN1N=CC=C1", "3-(1-Pyrazolyl)-alanine"),
        ("N[C@@H](CC1=CC=C(F)C=C1)C(O)=O", "4-fluorophenyl-L-alanine"),
        ("CSCCCCC(N(O)O)C(O)=O", "N,N-dihydroxydihomomethionine"),
        ("P(OC(C(N)C(O)=O)C)(O)(O)=O", "O-phosphonothreonine"),
        ("NC(Cc1c[nH]c2cccc(O)c12)C(O)=O", "4-hydroxytryptophan"),
        ("CSCCCCCC(NO)C(O)=O", "N-hydroxytrihomomethionine"),
        ("O(C[C@@H](N)C(O)=O)CC1=CC=CC=C1", "O-Benzyl-D-serine"),
        ("O(C1=CC=C(CCNC(=O)CCC(N)C(O)=O)C=C1)CC=2C=C(OC2)CN", "APMF-Glu"),
        ("CN[C@@H](Cc1ccccc1)C(O)=O", "N-methyl-L-phenylalanine"),
        ("C[C@H](NCCC(O)=O)C(O)=O", "(S)-beta-alanopine"),
        ("IC=1C=CC(C[C@H](C(O)=O)N)=CC1", "4-iodo-D-phenylalanine"),
        # Examples that were previously misclassified:
        ("O=C(N[C@@H](CC=1NC=NC1)C(O)=O)[C@@H](NC(=O)[C@@H](N)[C@H](O)C)CC2=CC=CC=C2", "Thr-Phe-His (peptide)"),
        ("Nc1cccc(c1)C(O)=O", "3-aminobenzoic acid"),
    ]
    for s, name in test_cases:
        result, reason = is_amino_acid(s)
        print(f"{name}: {result} -> {reason}")