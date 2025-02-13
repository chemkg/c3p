"""
Classifies: CHEBI:15841 polypeptide
"""
"""
Classifies: polypeptide – A peptide containing ten or more amino acid residues.
This improved version counts peptide bonds by ensuring the detected C(=O)N bond
belongs to a peptide backbone. It checks that:
  • The carbonyl carbon (C) has exactly three heavy neighbors: a double‐bonded oxygen,
    a nitrogen (the partner in the amide bond), and an “alpha” carbon (carbon attached
    to a side chain).
  • The amide nitrogen (N) has at least one other heavy neighbor that is a carbon (its
    associated alpha carbon).
If a valid peptide bond is found, then residue_count = (number of bonds) + 1.
If residue_count >= 10, the molecule is classified as a polypeptide.
Note: This is a heuristic and may still miss cases or mis‐classify very modified peptides
or other amide–rich compounds.
"""

from rdkit import Chem

def is_polypeptide(smiles: str):
    """
    Determines if a molecule is a polypeptide (a peptide with 10 or more amino acid residues)
    based on its SMILES string using an improved heuristic that verifies the bonding environment
    of each candidate peptide bond.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is judged to be a polypeptide, False otherwise.
        str: A message with the residue and peptide bond count details.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # SMARTS for an amide bond (C(=O)N); note that further filtering is done via environment tests.
    peptide_bond_pattern = Chem.MolFromSmarts("C(=O)N")
    if peptide_bond_pattern is None:
        return False, "Error creating SMARTS pattern"

    # Find all substructure matches for the amide bond.
    matches = mol.GetSubstructMatches(peptide_bond_pattern)
    valid_bond_count = 0

    for match in matches:
        c_idx, n_idx = match  # c_idx: carbonyl carbon, n_idx: amide nitrogen.
        c_atom = mol.GetAtomWithIdx(c_idx)
        n_atom = mol.GetAtomWithIdx(n_idx)

        # Check the environment of the carbonyl carbon.
        # In a peptide backbone the carbonyl carbon should be sp2 and have 3 heavy neighbors.
        neighbors_c = list(c_atom.GetNeighbors())
        if len(neighbors_c) != 3:
            continue  # Not typical peptide-backbone carbon.
        oxygen_found = False
        alphaC_found = False
        for neighbor in neighbors_c:
            bond = mol.GetBondBetweenAtoms(c_idx, neighbor.GetIdx())
            if neighbor.GetAtomicNum() == 8:
                # Look for a double bond to oxygen.
                if bond is not None and bond.GetBondTypeAsDouble() == 2.0:
                    oxygen_found = True
            elif neighbor.GetIdx() == n_idx:
                # This neighbor is the amide nitrogen (will inspect later).
                continue
            elif neighbor.GetAtomicNum() == 6:
                # An alpha carbon is expected.
                alphaC_found = True
        if not (oxygen_found and alphaC_found):
            continue

        # Check the environment of the amide nitrogen.
        # It should be bonded to the carbonyl carbon (already in match) and at least one alpha carbon.
        neighbors_n = list(n_atom.GetNeighbors())
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
        return True, f"Contains {residue_count} amino acid residues (found {valid_bond_count} peptide bonds)"
    else:
        return False, f"Only {residue_count} amino acid residues detected (found {valid_bond_count} peptide bonds)"

# Example usage for testing:
if __name__ == "__main__":
    # Test SMILES strings (use one of the provided examples):
    test_smiles = "O=C(N[C@H](C(=O)N[C@H](C(=O)O)CCCN(O)C(=O)C)CCCN(O)C(=O)C)[C@H]1N(C(=O)[C@H](NC(=O)[C@@H](NC(=O)[C@H]([C@H](O)[C@@H](NC(=O)CCCCCCCCCCC)CCCN=C(N)N)C)[C@H](O)C(=O)O)CO)CCC1"
    result, msg = is_polypeptide(test_smiles)
    print(result, msg)