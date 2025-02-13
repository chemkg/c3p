"""
Classifies: CHEBI:15705 L-alpha-amino acid
"""
"""
Classifies: L-alpha-amino acid
Definition: Any alpha-amino acid having L-configuration at the alpha-carbon.
The classification requires that the free (nonpeptidic) alpha-amino acid motif,
    N-[C@H](*)C(=O)[O;H1,O-]   or   N-[C@@H](*)C(=O)[O;H1,O-]
is present exactly once, with the amino nitrogen free (i.e., not amidated) and
the carboxyl group appearing as a free acid. The CIP code on the chiral alpha-carbon
must be 'S' (representing an L-amino acid).
"""

from rdkit import Chem

def is_L_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is an L-alpha amino acid from its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        (bool, str): True and a success reason if the molecule is a free L-alpha amino acid;
                     otherwise False and the reason for classification failure.
                     In ambiguous cases the function may return (None, None).
    """
    # Parse the input SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Ensure that stereochemistry and CIP codes are assigned properly
    Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
    
    # Define SMARTS patterns for the free alpha-amino acid motif.
    # We look for a nitrogen (as NH2 or NH3+), attached to a chiral carbon (either [C@H] or [C@@H])
    # which in turn is bound to a carboxyl group that is free (not engaged in peptide bonds).
    pattern1 = Chem.MolFromSmarts("[$(NH2),$(NH3+)]-[C@H](-*)(C(=O)[O;H1,O-])")
    pattern2 = Chem.MolFromSmarts("[$(NH2),$(NH3+)]-[C@@H](-*)(C(=O)[O;H1,O-])")
    
    matches1 = mol.GetSubstructMatches(pattern1)
    matches2 = mol.GetSubstructMatches(pattern2)
    all_matches = matches1 + matches2
    
    if not all_matches:
        return False, "No free alpha-amino acid motif found (expected pattern: N-[C@H](*)C(=O)[O;H1,O-])"
    
    # If more than one match is found the molecule might be a peptide (or polyamino acid)
    if len(all_matches) > 1:
        return False, f"Found {len(all_matches)} alpha-amino acid motifs; molecule likely represents a peptide"
    
    # Take the only match found
    match = all_matches[0]
    # According to our SMARTS, the indices are:
    #  index 0: amino nitrogen
    #  index 1: the chiral alpha–carbon
    #  index 2: the carboxyl carbon of the free acid group
    amino_idx = match[0]
    alpha_idx = match[1]
    carboxyl_idx = match[2]
    
    # Check that the amino nitrogen appears free (i.e., not amidated)
    amino_atom = mol.GetAtomWithIdx(amino_idx)
    # Count explicit + implicit hydrogens on the amino nitrogen
    total_H = amino_atom.GetTotalNumHs()
    if total_H < 2:
        return False, f"Amino nitrogen has too few hydrogens ({total_H}); may be amidated"
    
    # Additional check: ensure that aside from its bond to the alpha-carbon, the nitrogen does not bond to a carbonyl carbon.
    for nbr in amino_atom.GetNeighbors():
        if nbr.GetIdx() == alpha_idx:
            continue
        if nbr.GetAtomicNum() == 6:  # carbon
            for bond in nbr.GetBonds():
                # Look for a double bond from this carbon to oxygen
                if bond.GetBondTypeAsDouble() == 2.0:
                    other = bond.GetOtherAtom(nbr)
                    if other.GetAtomicNum() == 8:
                        return False, "Amino nitrogen appears amidated (neighboring carbonyl group detected)"
    
    # Check the carboxyl group is a free acid.
    carboxyl_atom = mol.GetAtomWithIdx(carboxyl_idx)
    o_neighbors = [nbr for nbr in carboxyl_atom.GetNeighbors() if nbr.GetAtomicNum() == 8]
    if len(o_neighbors) != 2:
        return False, f"Carboxyl carbon does not have exactly 2 oxygen neighbors (found {len(o_neighbors)})"
    
    # Ensure neither oxygen is bonded to a nitrogen (which may indicate amidation)
    for o in o_neighbors:
        for nbr in o.GetNeighbors():
            if nbr.GetAtomicNum() == 7 and nbr.GetIdx() != carboxyl_idx:
                return False, "Carboxyl group appears amidated (oxygen bonded to nitrogen)"
    
    # Optional sanity check: if the molecule contains multiple amide bonds beyond the free acid, it may be a peptide.
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    if len(mol.GetSubstructMatches(amide_pattern)) > 1:
        return False, "Multiple amide bonds detected – molecule likely represents a peptide"
    
    # Check the CIP configuration of the alpha-carbon.
    alpha_atom = mol.GetAtomWithIdx(alpha_idx)
    if not alpha_atom.HasProp('_CIPCode'):
        return False, "Alpha-carbon chiral configuration is not assigned"
    actual_config = alpha_atom.GetProp('_CIPCode')
    expected_config = 'S'  # For L-amino acids under our convention
    if actual_config != expected_config:
        return False, f"Alpha-carbon configuration is {actual_config}, expected {expected_config} for L-amino acid"
    
    return True, "Contains free L-alpha-amino acid motif with correct configuration"

# Example usage:
if __name__ == "__main__":
    # Test examples: feel free to add more from the provided list.
    test_smiles = [
        "CC(C)(CCO)SC[C@H](N)C(O)=O",  # felinine
        "CC(=O)[C@H](N)C(O)=O",  # L-2-amino-3-oxobutanoic acid
        "N[C@@H](CSCCB(O)O)C(O)=O",  # S-(2-boronoethyl)-L-cysteine
        "CC(C)[C@H](N)C(O)=O",  # L-valine
        "N[C@@H](CCCCNC(O)=O)C(O)=O"  # N(6)-carboxy-L-lysine
    ]
    for s in test_smiles:
        res, reason = is_L_alpha_amino_acid(s)
        print(s, "=>", res, "|", reason)