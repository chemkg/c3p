"""
Classifies: CHEBI:87659 dodecanoate ester
"""
#!/usr/bin/env python3
"""
Classifies: dodecanoate ester
Defined as: Any fatty acid ester in which the carboxylic acid component is lauric acid.
Lauric acid (dodecanoic acid) is CH3(CH2)10COOH; when esterified the acyl moiety is
CH3(CH2)10C(=O)O–. Hence, we search for an ester group that exactly encodes a laurate residue.

Algorithm:
  1. Parse the SMILES.
  2. Use a SMARTS that “hard‐codes” the laurate acyl group using proper repetition:
       "[CH3]-([CH2]){10}-C(=O)O[C]"
     This means a terminal CH3, followed by exactly 10 CH2’s, then a carbonyl (C(=O)) 
     and an oxygen that is attached to a carbon.
  3. For each substructure match found, check that the first atom (the CH3) is terminal; 
     that is, it should have exactly one bonded carbon that is part of the match.
  4. If any match passes this test, report the molecule as containing a dodecanoate ester.
"""

from rdkit import Chem

def is_dodecanoate_ester(smiles: str):
    """
    Determines if a molecule (given by its SMILES string) is a dodecanoate ester.
    We require that the molecule contains an ester group in which the acyl part
    is exactly lauric acid, CH3(CH2)10C(=O)O – where the terminal CH3 is isolated.

    Args:
      smiles (str): A SMILES string.

    Returns:
      (bool, str): A tuple with classification and a reason.
    """
    # Parse the SMILES string to a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."

    # Build the SMARTS query for a laurate ester.
    # Use parentheses for the repeated [CH2] group.
    laurate_smarts = "[CH3]-([CH2]){10}-C(=O)O[C]"
    query = Chem.MolFromSmarts(laurate_smarts)
    if query is None:
        return False, "Error creating SMARTS for laurate ester."

    # Find all substructure matches.
    matches = mol.GetSubstructMatches(query)
    if not matches:
        return False, "Molecule does not contain a dodecanoate (laurate) ester moiety."

    # For each match, confirm that the first matched atom (the terminal CH3) is isolated.
    for match in matches:
        # The first atom in our SMARTS is the terminal CH3.
        first_idx = match[0]
        first_atom = mol.GetAtomWithIdx(first_idx)
        # Check that this CH3 does not have any additional carbon neighbors outside the match.
        # A terminal methyl should only have one carbon neighbor.
        carbon_neighbors = [nbr for nbr in first_atom.GetNeighbors() if nbr.GetAtomicNum() == 6]
        # If there is not exactly one carbon neighbor, then this CH3 may be extended.
        if len(carbon_neighbors) != 1:
            continue
        # Also ensure that the single carbon neighbor is part of the match.
        if carbon_neighbors[0].GetIdx() not in match:
            continue
        # If we pass all conditions, we have found a dodecanoate ester.
        return True, "Molecule contains a dodecanoate (laurate) ester functionality."

    # If no valid match is found, the acyl chain may be extended.
    return False, ("Molecule contains a fragment matching laurate ester connectivity, "
                   "but the acyl chain is extended (the terminal CH3 is not isolated).")

# Example usage:
if __name__ == "__main__":
    # Some test SMILES strings.
    test_smiles_list = [
        "CCCCCCCCCCCC(=O)OCC",  # simple dodecanoate ester (laurate ester)
        "O(CCCCCC(C)C)C(=O)CCCCCCCCCCCC",  # another example; may be more complex
        "invalid_smiles"
    ]
    for s in test_smiles_list:
        result, reason = is_dodecanoate_ester(s)
        print(f"SMILES: {s}\nResult: {result}\nReason: {reason}\n")