"""
Classifies: CHEBI:52221 isothiocyanate
"""
"""
Classifies: Isothiocyanate (organosulfur compound with the general formula R-N=C=S)

Improved version:
- Uses an explicit SMARTS pattern to capture an sp2 nitrogen double‐bonded to an sp2 carbon,
  which in turn is double–bonded to a terminal sulfur.
- Checks that each atom in the group has the expected connectivity.
- Checks that the R substituent at nitrogen is not immediately part of an amide (i.e. not directly attached
  to a carbonyl C=O) to filter out many false positives.
"""

from rdkit import Chem

def is_isothiocyanate(smiles: str):
    """
    Determines if a molecule is a (simple) isothiocyanate based on its SMILES string.
    A simple isothiocyanate is defined as an organosulfur compound of the form R-N=C=S,
    where the N, C, and S atoms have the expected bonding pattern and the substituent R on N is not
    immediately part of an amide carbonyl motif.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule contains a valid isothiocyanate group, False otherwise.
        str: Reason for the classification.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Use an explicit SMARTS pattern that requires:
    # - A nitrogen (sp2) double-bonded to a carbon (sp2)
    # - That carbon is double-bonded to a sulfur that is not further substituted.
    # Note: We do not force the R group at nitrogen in the SMARTS since it can be any atom.
    iso_smarts = "[NX2]=[CX2]=[SX1]"
    iso_group = Chem.MolFromSmarts(iso_smarts)
    
    if iso_group is None:
        return False, "Error in generating isothiocyanate SMARTS"
    
    # If no isothiocyanate substructure is found, reject.
    if not mol.HasSubstructMatch(iso_group):
        return False, "No isothiocyanate functional group (N=C=S) found"
    
    # Get all substructure matches.
    matches = mol.GetSubstructMatches(iso_group)
    valid_match_found = False
    for match in matches:
        # match is a tuple of atom indices corresponding to N, C, S (in the order as defined in SMARTS).
        n_idx, c_idx, s_idx = match
        n_atom = mol.GetAtomWithIdx(n_idx)
        c_atom = mol.GetAtomWithIdx(c_idx)
        s_atom = mol.GetAtomWithIdx(s_idx)
        
        # Check connectivity:
        # For a valid simple isothiocyanate we expect:
        # - Nitrogen: bonded to the matched carbon and one additional substituent (the R group)
        # - Carbon: bonded only to the N and S in this group.
        # - Sulfur: terminal (bonded only to the carbon).
        if n_atom.GetDegree() != 2:
            continue  # N must connect only to the isothiocyanate carbon and one more atom (R)
        if c_atom.GetDegree() != 2:
            continue  # C should be only connected to N and S.
        if s_atom.GetDegree() != 1:
            continue  # S must be terminal.
        
        # Also check that formal charges are 0 for all three.
        if n_atom.GetFormalCharge() != 0 or c_atom.GetFormalCharge() != 0 or s_atom.GetFormalCharge() != 0:
            continue
        
        # Retrieve the substituent (R group) attached to the nitrogen (the neighbor that is not the matched carbon).
        n_neighbors = [nbr for nbr in n_atom.GetNeighbors() if nbr.GetIdx() != c_idx]
        if not n_neighbors:
            continue  # There must be a substituent.
        r_atom = n_neighbors[0]  # Considering only one substituent.
        
        # As an extra filtering step, check that the R substituent is not immediately attached to a carbonyl,
        # which is typical of amide or related motifs known to cause false positives.
        r_bad = False
        for nbr in r_atom.GetNeighbors():
            # Look for a bond from r_atom to a carbon that is double-bonded to oxygen:
            if nbr.GetAtomicNum() == 6:
                for bond in nbr.GetBonds():
                    # Look for a double bond to oxygen.
                    if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                        other = bond.GetOtherAtom(nbr)
                        if other.GetAtomicNum() == 8:
                            r_bad = True
                            break
            if r_bad:
                break
        if r_bad:
            continue  # Skip this match if the substituent appears to be part of an amide-like motif.
        
        # If one valid match is found, we consider the molecule as an isothiocyanate.
        valid_match_found = True
        break

    if valid_match_found:
        return True, "Contains isothiocyanate functional group (R-N=C=S) with appropriate substitution"
    else:
        return False, "Isothiocyanate group found does not meet criteria for a simple R-N=C=S motif"
        
# You can test the function with some examples.
if __name__ == "__main__":
    test_smiles = [
        "S(CCCCCCCCN=C=S)C",  # 8-(methylthio)octylisothiocyanate, expected True
        "S(=O)(CCCCCCCN=C=S)C"  # 1-isothiocyanato-7-(methylsulfinyl)heptane, known FP; expected False
    ]
    for smi in test_smiles:
        result, reason = is_isothiocyanate(smi)
        print(f"SMILES: {smi}\nResult: {result}\nReason: {reason}\n")