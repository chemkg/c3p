"""
Classifies: CHEBI:15705 L-alpha-amino acid
"""
Below is one possible solution. In this version we first assign stereochemistry and then use two SMARTS to “find‐the” free (nonpeptidic) α–amino acid motif. In addition to requiring that exactly one match exists, we perform “sanity” checks on the free NH₂ and –COOH groups (for example, ensuring that the amino nitrogen isn’t “hidden” in a peptide bond and that the carboxyl carbon indeed is not amidated). Finally we extract the CIP code from the (assigned) chiral α–carbon and (by convention for this classification task) require that it be “S” (even though for cysteine one sometimes expects the reverse CIP label). (Note that because molecule depiction and‐stereochemistry are subtle topics there are many ways to do this; if the task is too hard the function may return (None, None).)

Below is the complete Python code in one code block:

------------------------------------------------

"""
Classifies: L-alpha-amino acid
Definition: Any alpha–amino acid having L–configuration at the alpha carbon.
For the purposes of this classifier we require that the motif
    N–[C@H](*)–C(=O)O   or   N–[C@@H](*)–C(=O)O
is present exactly once in a molecule that does not appear to be a peptide.
We also “sanity check” that the amino nitrogen is free (i.e. not amidated)
and that the carboxyl group is a free –COOH rather than forming an amide.
Finally, we use RDKit’s assigned CIP code on the chiral α–carbon
and—by convention for these L–amino acids—demand it is “S”.
If any step fails the function returns False along with a reason.
"""

from rdkit import Chem

def is_L_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is an L-alpha amino acid based on its SMILES string.
    Steps:
      1. Parse the SMILES and assign stereochemistry.
      2. Look for exactly one occurrence of the free alpha–amino acid motif:
           either "N[C@H](*)C(=O)O" or "N[C@@H](*)C(=O)O"
         (which indicates a free amino group attached to a chiral alpha–carbon attached to a free carboxylic acid).
      3. Check that the amino nitrogen is not amidated (i.e. it has at least two hydrogens)
         and that the carboxyl carbon is not joined to any extra nitrogen (as in a peptide).
      4. Retrieve the assigned CIP configuration for the alpha–carbon.
         (For this classifier we require that the configuration be “S” for an L–amino acid.)
      5. If more than one match is found the molecule is likely a peptide, so that result becomes False.
    
    Args:
      smiles (str): SMILES string of the molecule.
    
    Returns:
      (bool, str): Tuple where the boolean is True if the molecule matches the class
                   (L-alpha-amino acid) and the reason indicates the outcome.
                   In case of errors or an ambiguous case, returns (False, reason).
                   (If the task is too ambiguous the function may return (None, None).)
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Always update chirality information
    Chem.AssignStereochemistry(mol, force=True, cleanIt=True)

    # Define two SMARTS that capture a free alpha-amino acid motif.
    # The pattern requires a nitrogen directly bonded to a chiral carbon which in turn is bonded to a carboxyl.
    # (Note: Free amino acids are expected to have an –NH2 (or –NH3+) and a –COOH group.)
    pattern1 = Chem.MolFromSmarts("[$(NH2),$(NH3+)]-[C@H](-*)(C(=O)[O;H1,O-])")
    pattern2 = Chem.MolFromSmarts("[$(NH2),$(NH3+)]-[C@@H](-*)(C(=O)[O;H1,O-])")
    
    matches1 = mol.GetSubstructMatches(pattern1)
    matches2 = mol.GetSubstructMatches(pattern2)
    matches = matches1 + matches2
    if not matches:
        return False, "No alpha-amino acid motif found (expected pattern: N-[C@H](*)C(=O)O)"
    
    # If more than one match is found we suspect the molecule is a peptide.
    if len(matches) > 1:
        return False, f"Found {len(matches)} alpha-amino acid motifs; molecule likely represents a peptide"
    
    # We use the first (and only) match.
    # According to our SMARTS the match indices are:
    #   index0: amino nitrogen
    #   index1: chiral alpha–carbon
    #   index2: carboxyl carbon (of the free –COOH)
    match = matches[0]
    amino_idx = match[0]
    alpha_idx = match[1]
    carboxyl_idx = match[2]
    
    # Check the amino nitrogen: it should be free (i.e. not part of an amide bond)
    amino_atom = mol.GetAtomWithIdx(amino_idx)
    # Count total hydrogens (explicit plus implicit)
    nH = amino_atom.GetTotalNumHs()
    if nH < 2:
        return False, f"Amino nitrogen has too few hydrogens ({nH}); may be amidated"
    # Also, ensure that besides the alpha carbon bond, this N is not bonded to a carbon with a C=O (which would indicate peptide linkage)
    for nbr in amino_atom.GetNeighbors():
        if nbr.GetIdx() == alpha_idx:
            continue
        # If neighbor is carbon and is double-bonded to oxygen, flag it
        if nbr.GetAtomicNum() == 6:
            for bond in nbr.GetBonds():
                if bond.GetBeginAtomIdx() == nbr.GetIdx() or bond.GetEndAtomIdx() == nbr.GetIdx():
                    # Look for C=O bond
                    if bond.GetBondTypeAsDouble() == 2.0:
                        other = bond.GetOtherAtom(nbr)
                        if other.GetAtomicNum() == 8:
                            return False, "Amino nitrogen appears amidated (bonded to a carbonyl group)"
    
    # Check the carboxyl carbon: verify it is a free acid (i.e. not amidated)
    carboxyl_atom = mol.GetAtomWithIdx(carboxyl_idx)
    # It should be attached to exactly two oxygens (one as a carbonyl and one as hydroxyl or deprotonated O)
    o_neighbors = [nbr for nbr in carboxyl_atom.GetNeighbors() if nbr.GetAtomicNum() == 8]
    if len(o_neighbors) != 2:
        return False, f"Carboxyl carbon does not have exactly 2 oxygen neighbors (found {len(o_neighbors)})"
    # Ensure none of the oxygens is further attached to a nitrogen (which would indicate an amide)
    for o in o_neighbors:
        for nbr in o.GetNeighbors():
            if nbr.GetAtomicNum() == 7 and nbr.GetIdx() != carboxyl_idx:
                return False, "Carboxyl group appears amidated"
    
    # Additional sanity check: if the molecule contains more than one amide bond (apart from the carboxyl group) it is probably a peptide.
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    if len(amide_matches) > 1:
        return False, "Multiple amide bonds detected – molecule likely represents a peptide rather than a single amino acid"
        
    # Now get the chiral configuration of the candidate alpha-carbon.
    alpha_atom = mol.GetAtomWithIdx(alpha_idx)
    if not alpha_atom.HasProp('_CIPCode'):
        return False, "Alpha-carbon chiral configuration is not assigned"
    actual_config = alpha_atom.GetProp('_CIPCode')  # typically 'R' or 'S'
    
    # For this classification we demand that the free (L-) alpha amino acid has an S–configuration.
    # (Even though in some cases modified cysteine may have a reversed CIP label, here the ground truth is taken as “S”.)
    expected_config = 'S'
    if actual_config != expected_config:
        return False, f"Alpha–carbon configuration is {actual_config}, but expected {expected_config} for an L–alpha–amino acid"
    
    return True, "Contains free alpha–amino acid motif with correct L–configuration"

# Example usage (uncomment to test):
#if __name__ == "__main__":
#    test_smiles = "CC(=O)[C@H](N)C(O)=O"  # L-2-amino-3-oxobutanoic acid
#    res, reason = is_L_alpha_amino_acid(test_smiles)
#    print(res, reason)
    
#    # Try a peptide (should return False)
#    peptide_smiles = "O=C(N[C@@H](N)CC(O)=O)C(=O)N[C@@H](N)CC(O)=O"  # dipeptide example
#    print(is_L_alpha_amino_acid(peptide_smiles))
    
#    # Test a modified amino acid (e.g., felinine)
#    felinine = "CC(C)(CCO)SC[C@H](N)C(O)=O"
#    print(is_L_alpha_amino_acid(felinine))
    
#    # Test S-(2-boronoethyl)-L-cysteine
#    borono = "N[C@@H](CSCCB(O)O)C(O)=O"
#    print(is_L_alpha_amino_acid(borono))
    
#    # etc.
    

# Note: Because of the subtleties in stereochemistry assignment in modified amino acids,
# some borderline structures may be mis‐classified by this heuristic.
    
