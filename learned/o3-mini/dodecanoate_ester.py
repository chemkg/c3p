"""
Classifies: CHEBI:87659 dodecanoate ester
"""
#!/usr/bin/env python3
"""
Classifies: dodecanoate ester
Defined as: Any fatty acid ester in which the carboxylic acid component is lauric acid.
Lauric acid (dodecanoic acid) is CH3(CH2)10COOH; when esterified the acyl moiety is
CH3(CH2)10C(=O)O–. This program searches for a SMARTS match that exactly encodes a laurate residue.
It then verifies that the terminal CH3 is not “extended” (i.e. that it is not a part of a larger chain).
"""

from rdkit import Chem

def is_dodecanoate_ester(smiles: str):
    """
    Determines if a molecule (given by its SMILES string) is a dodecanoate ester.
    We require that the molecule contains an ester group in which the acyl (fatty acid) part
    is exactly lauric acid, CH3(CH2)10C(=O)O– (i.e. 12 carbons, where the first atom is a terminal CH3).
    
    The algorithm is:
      1. Parse the SMILES.
      2. Use a SMARTS that “hard‐codes” the laurate acyl group:
           "[CH3]-[CH2]{10}-C(=O)O[C]"
         This means a CH3, followed by exactly 10 CH2’s then a carbonyl (making up 12 carbons),
         and the ester oxygen must be attached to a carbon.
      3. For each substructure match found, check that the first atom of the match (the CH3)
         is terminal, i.e. none of its other carbon neighbours lies outside the match:
             That is, the CH3 should have exactly 1 carbon neighbor and that neighbor is also in the match.
      4. If any match passes this test, report the molecule as containing a dodecanoate ester.
    
    Args:
      smiles (str): A SMILES string.
    
    Returns:
      (bool, str): A tuple with True/False and a reason.
    """
    # Parse the SMILES string to a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."

    # Build the SMARTS query for a laurate ester.
    # The pattern: terminal CH3, then exactly 10 CH2 groups, then a carbonyl (C(=O)) linked to O attached to a carbon.
    laurate_smarts = "[CH3]-[CH2]{10}-C(=O)O[C]"
    query = Chem.MolFromSmarts(laurate_smarts)
    if query is None:
        return False, "Error creating SMARTS for laurate ester."

    # Find all substructure matches.
    matches = mol.GetSubstructMatches(query)
    if not matches:
        return False, "Molecule does not contain a dodecanoate (laurate) ester moiety."

    # For each match, confirm that the first matched atom (the CH3) is terminal.
    for match in matches:
        # match is a tuple of atom indices. The order is the same as in the SMARTS.
        # The first atom should be the terminal CH3.
        first_idx = match[0]
        first_atom = mol.GetAtomWithIdx(first_idx)
        # Check that this atom is a methyl: degree for carbon should be exactly 1 (the chain connection)
        # Count bonded carbon neighbors that are not in the match.
        external_carbon = False
        for nbr in first_atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6 and nbr.GetIdx() not in match:
                external_carbon = True
                break
        if external_carbon:
            continue  # This match is part of a longer chain than laurate.
        # If a valid match is found then we classify as a dodecanoate ester.
        return True, "Molecule contains a dodecanoate (laurate) ester functionality."

    return False, ("Molecule contains a fragment matching laurate ester connectivity, "
                   "but the acyl chain is extended (the terminal CH3 is not isolated).")

# Example usage:
if __name__ == "__main__":
    # Some short tests – you can add more SMILES examples for testing.
    test_smiles_list = [
        "CCCCCCCCCCCC(=O)OCC",  # simple dodecanoate ester (laurate ester)
        "O(CCCCCC(C)C)C(=O)CCCCCCCCCCCC",  # example (might be more complex)
        "invalid_smiles"
    ]
    for s in test_smiles_list:
        result, reason = is_dodecanoate_ester(s)
        print(f"SMILES: {s}\nResult: {result}\nReason: {reason}\n")