"""
Classifies: CHEBI:15705 L-alpha-amino acid
"""
"""
Classifies: L-alpha-amino acid
Definition: Any alpha-amino acid having L-configuration at the alpha-carbon.
A free (nonpeptidic) L-alpha-amino acid should have exactly one motif of the form:
    [NH2,NH3+]-[C@H or C@@H]([*])C(=O)[O or O-]
with the amino nitrogen free (not amidated) and the carboxyl group as a free acid.
In addition, the CIP code for the chiral (alpha) carbon should be 'S' (L configuration).
"""

from rdkit import Chem

def is_L_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is a free L-alpha amino acid from its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        (bool, str): (True, reason) if the molecule is a free L-alpha amino acid,
                     (False, reason) if not, or (None, None) in ambiguous cases.
    """
    # Parse the input SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Ensure stereochemistry is assigned so that CIP codes become available.
    Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
    
    # Define SMARTS patterns for the free alpha-amino acid motif.
    # We allow two forms for the carboxyl group: protonated (-C(=O)O) and deprotonated (-C(=O)[O-]),
    # and two notations for the chiral alpha carbon: [C@H] or [C@@H].
    # The amino group is specified as either NH2 or NH3+.
    patterns = []
    smarts_list = [
        "[NH2,NH3+]-[C@H]([*])C(=O)[O-]",  # protonated N, chiral as [C@H], carboxyl deprotonated
        "[NH2,NH3+]-[C@H]([*])C(=O)O",     # protonated N, chiral as [C@H], carboxyl protonated
        "[NH2,NH3+]-[C@@H]([*])C(=O)[O-]", # protonated N, chiral as [C@@H], carboxyl deprotonated
        "[NH2,NH3+]-[C@@H]([*])C(=O)O"     # protonated N, chiral as [C@@H], carboxyl protonated
    ]
    for smarts in smarts_list:
        patt = Chem.MolFromSmarts(smarts)
        if patt is None:
            # If any pattern fails to compile, we cannot proceed reliably.
            return None, None
        patterns.append(patt)
    
    # Find all matches for any of the defined patterns.
    all_matches = []
    for patt in patterns:
        matches = mol.GetSubstructMatches(patt, useChirality=True)
        for m in matches:
            if m not in all_matches:
                all_matches.append(m)
    
    if not all_matches:
        return False, "No free alpha-amino acid motif found; expected pattern: [NH2,NH3+]-[C@H/C@@H]([*])C(=O)[O or O-]"
    
    # If more than one motif is found, the molecule is likely a peptide.
    if len(all_matches) > 1:
        return False, f"Found {len(all_matches)} alpha-amino acid motifs; molecule likely represents a peptide"
    
    # Using the single match found, extract the indices.
    # Based on our SMARTS, the mapping is:
    #  index 0: amino nitrogen
    #  index 1: the chiral (alpha) carbon
    #  index 2: the carboxyl carbon
    match = all_matches[0]
    amino_idx = match[0]
    alpha_idx = match[1]
    carboxyl_idx = match[2]
    
    # Check that the amino nitrogen appears "free" (i.e., it should have at least 2 hydrogens).
    amino_atom = mol.GetAtomWithIdx(amino_idx)
    total_H = amino_atom.GetTotalNumHs()
    if total_H < 2:
        return False, f"Amino nitrogen has too few hydrogens ({total_H}); may be amidated"
    
    # Ensure the carboxyl carbon is correctly connected:
    # It should have one double bond oxygen and one single bond oxygen (hydroxyl or deprotonated oxygen).
    carboxyl_atom = mol.GetAtomWithIdx(carboxyl_idx)
    oxygens = [nbr for nbr in carboxyl_atom.GetNeighbors() if nbr.GetAtomicNum() == 8]
    if len(oxygens) != 2:
        return False, f"Carboxyl carbon does not have exactly 2 oxygen neighbors (found {len(oxygens)})"
    double_bond_found = False
    single_bond_found = False
    for o in oxygens:
        bond = mol.GetBondBetweenAtoms(carboxyl_idx, o.GetIdx())
        if bond is None:
            continue
        # Check if the bond is a double bond (representing C=O)
        if bond.GetBondTypeAsDouble() == 2.0:
            double_bond_found = True
        else:
            single_bond_found = True
    if not double_bond_found or not single_bond_found:
        return False, "Carboxyl group does not appear to be a free acid (expected one double and one single bond to O)"
    
    # Check for extra amide bonds that might indicate a peptide.
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    if len(amide_matches) > 1:
        return False, "Multiple amide bonds detected – molecule likely represents a peptide"
    
    # Check that the chiral (alpha) carbon has been assigned a CIP code.
    alpha_atom = mol.GetAtomWithIdx(alpha_idx)
    if not alpha_atom.HasProp('_CIPCode'):
        return False, "Alpha-carbon stereochemistry (CIP code) is not assigned"
    actual_config = alpha_atom.GetProp('_CIPCode')
    expected_config = 'S'
    if actual_config != expected_config:
        return False, f"Alpha-carbon configuration is {actual_config}, expected {expected_config} for L-amino acid"
    
    return True, "Contains free L-alpha-amino acid motif with correct configuration"

# Example usage (testing with several provided SMILES strings):
if __name__ == "__main__":
    test_smiles = [
        "CC(C)(CCO)SC[C@H](N)C(O)=O",  # felinine
        "CC(=O)[C@H](N)C(O)=O",        # L-2-amino-3-oxobutanoic acid
        "N[C@@H](CSCCB(O)O)C(O)=O",     # S-(2-boronoethyl)-L-cysteine
        "N[C@@H](CCCCNC(O)=O)C(O)=O",   # N(6)-carboxy-L-lysine
        "CC(C)[C@H](N)C(O)=O"          # L-valine
    ]
    for s in test_smiles:
        res, reason = is_L_alpha_amino_acid(s)
        print(s, "=>", res, "|", reason)