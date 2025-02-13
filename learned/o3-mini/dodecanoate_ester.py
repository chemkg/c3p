"""
Classifies: CHEBI:87659 dodecanoate ester
"""
#!/usr/bin/env python3
"""
Classifies: dodecanoate ester
Defined as: Any fatty acid ester in which the carboxylic acid component is lauric acid.
Lauric acid (dodecanoic acid) is CH3(CH2)10COOH. When esterified the acyl part becomes
CH3(CH2)10C(=O)O– (attached to an alcohol). The idea is to recognize an acyl fragment
with exactly 12 carbons (one CH3 + 10 CH2 + a carbonyl C) that is linked via an ester bond.
This improved program first looks for a SMARTS match for an ester fragment that contains
the laurate acyl group and then verifies that the first (terminal) carbon is not
“extended” (i.e. not part of a longer chain).
"""

from rdkit import Chem

def is_dodecanoate_ester(smiles: str):
    """
    Determines if a molecule (given by its SMILES string) is a dodecanoate ester.
    The molecule must contain an ester linkage in which the acyl (fatty acid) part is
    exactly lauric acid (dodecanoic acid): CH3(CH2)10C(=O)O–.
    
    The algorithm is as follows:
      1. Parse the SMILES string.
      2. Search for a substructure matching a laurate ester fragment.
         We build a SMARTS that explicitly encodes the acyl chain:
         We want the acyl group to be exactly 12 carbons in a row ending in a carboxylate part.
         In lauric acid ester the fragment appears as:
           CH3-(CH2)10-C(=O)O–  (attached to some R group).
         We build the SMARTS by first writing out 11 consecutive aliphatic carbons (which,
         when followed by one more “C(=O)” gives 12 carbons total) and then appending the carbonyl and ester oxygen.
         Finally, we require that the oxygen is bonded to a carbon (and not, say, to a non‐carbon atom).
         The resulting SMARTS is: "CCCCCCCCCCCC(=O)O[C]"
         (Note: When concatenated, the pattern “CCCCCCCCCCC” (11 C’s) plus “C(=O)O[C]” gives a total
         of 12 carbons in the acyl fragment.)
      3. For each substructure match found, check that the very first matched carbon is terminal
         (i.e. it does not have a neighboring carbon outside the match which would indicate the chain continues).
      4. If any match passes that “end‐group” test, return True with an explanation.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a dodecanoate ester, False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # Define a SMARTS pattern for a laurate ester fragment.
    # Explanation:
    #   "CCCCCCCCCCC" is 11 consecutive aliphatic carbons.
    #   When we append "C(=O)O[C]" the extra C becomes the carbonyl carbon, making the acyl group exactly 12 carbons.
    #   The "[C]" after the ester oxygen ensures that the oxygen is linked to a carbon (the alcohol part).
    laurate_smarts = "CCCCCCCCCCC" + "C(=O)O[C]"
    query = Chem.MolFromSmarts(laurate_smarts)
    if query is None:
        return False, "Error creating SMARTS for laurate ester."
    
    # Get all substructure matches.
    matches = mol.GetSubstructMatches(query)
    if not matches:
        return False, "Molecule does not contain a dodecanoate (laurate) ester moiety."
    
    # For each match, ensure the acyl chain is not part of a longer chain.
    # The first atom of the match (match[0]) is expected to be the terminal CH3 group.
    # Check that none of its carbon neighbors lies outside of the matched substructure.
    for match in matches:
        # match is a tuple of atom indices corresponding to the SMARTS pattern.
        first_atom = mol.GetAtomWithIdx(match[0])
        valid = True
        for nbr in first_atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6 and nbr.GetIdx() not in match:
                # The terminal CH3 should not be connected to another carbon (i.e. it must be terminal).
                valid = False
                break
        if valid:
            return True, "Molecule contains a dodecanoate (laurate) ester functionality."
    
    return False, "Molecule contains a fragment matching laurate ester connectivity but it seems embedded in a longer chain."

# Example usage:
if __name__ == "__main__":
    # A test example: one known laurate ester fragment.
    test_smiles = "O(CCCCCC(C)C)C(=O)CCCCCCCCCCCC"  # This one is only for demo; real examples may be more complex.
    result, reason = is_dodecanoate_ester(test_smiles)
    print(f"Result: {result}\nReason: {reason}")