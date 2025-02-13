"""
Classifies: CHEBI:15841 polypeptide
"""
"""
Classifies: polypeptide – A peptide containing ten or more amino acid residues.
This version counts peptide bonds by verifying the bonding environment of each candidate amide bond,
and it carefully unpacks substructure matches to avoid "too many values to unpack" errors.
"""

from rdkit import Chem

def is_polypeptide(smiles: str):
    """
    Determines if a molecule is a polypeptide (a peptide with 10 or more amino acid residues)
    based on its SMILES string. This function uses an improved heuristic that tests the bonding 
    environment around each candidate peptide (amide) bond.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        (bool, str): A tuple with a boolean indicating whether the molecule is a polypeptide,
                     and a message with details on the residue and peptide bond counts.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Create SMARTS pattern for an amide bond (C(=O)N).
    peptide_bond_pattern = Chem.MolFromSmarts("C(=O)N")
    if peptide_bond_pattern is None:
        return False, "Error creating SMARTS pattern"
    
    # Find all substructure matches for the amide bond.
    matches = mol.GetSubstructMatches(peptide_bond_pattern)
    valid_bond_count = 0

    for match in matches:
        # Ensure match length is at least 2.
        if len(match) < 2:
            continue  # Skip invalid match tuple.
        # If tuple has more than 2 atoms, select the first two indices.
        if len(match) != 2:
            c_idx = match[0]
            n_idx = match[1]
        else:
            c_idx, n_idx = match

        try:
            c_atom = mol.GetAtomWithIdx(c_idx)
            n_atom = mol.GetAtomWithIdx(n_idx)
        except Exception:
            continue  # Skip if indices are not valid.

        # Check the environment of the carbonyl carbon:
        # - Should be sp2 and have exactly three heavy-atom neighbors.
        neighbors_c = [atom for atom in c_atom.GetNeighbors() if atom.GetAtomicNum() > 1]
        if len(neighbors_c) != 3:
            continue  # Not typical for a peptide-backbone carbon.
        
        oxygen_found = False
        alphaC_found = False
        for neighbor in neighbors_c:
            bond = mol.GetBondBetweenAtoms(c_idx, neighbor.GetIdx())
            # Look for a double bond to oxygen.
            if neighbor.GetAtomicNum() == 8:
                if bond is not None and bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                    oxygen_found = True
            elif neighbor.GetIdx() == n_idx:
                continue  # This neighbor is the amide nitrogen; will be checked separately.
            elif neighbor.GetAtomicNum() == 6:
                alphaC_found = True
        if not (oxygen_found and alphaC_found):
            continue

        # Check the environment of the amide nitrogen:
        # - Should be bonded to the carbonyl carbon and at least one other carbon (alpha carbon).
        neighbors_n = [atom for atom in n_atom.GetNeighbors() if atom.GetAtomicNum() > 1]
        if len(neighbors_n) < 2:
            continue
        alphaN_found = False
        for neighbor in neighbors_n:
            if neighbor.GetIdx() == c_idx:
                continue
            if neighbor.GetAtomicNum() == 6:
                alphaN_found = True
                break
        if not alphaN_found:
            continue

        # This match appears to be a valid peptide bond.
        valid_bond_count += 1

    # For a linear peptide the number of residues = (peptide bonds) + 1.
    residue_count = valid_bond_count + 1

    if residue_count >= 10:
        return True, f"Contains {residue_count} amino acid residues (found {valid_bond_count} peptide bonds)."
    else:
        return False, f"Only {residue_count} amino acid residues detected (found {valid_bond_count} peptide bonds)."

# Example usage for testing:
if __name__ == "__main__":
    # Test using a SMILES string that previously caused a "too many values to unpack" error.
    test_smiles = "O=C(N[C@H](C(=O)N[C@H](C(=O)O)CCCN(O)C(=O)C)CCCN(O)C(=O)C)[C@H]1N(C(=O)[C@H](NC(=O)[C@@H](NC(=O)[C@H]([C@H](O)[C@@H](NC(=O)CCCCCCCCCCC)CCCN=C(N)N)C)[C@H](O)C(=O)O)CO)CCC1"
    result, message = is_polypeptide(test_smiles)
    print(result, message)