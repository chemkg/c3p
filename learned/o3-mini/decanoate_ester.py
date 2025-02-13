"""
Classifies: CHEBI:87658 decanoate ester
"""
"""
Classifies: Decanoate ester

A decanoate ester is defined as a fatty acid ester resulting from the formal condensation 
of the carboxy group of decanoic (capric) acid with the hydroxy group of an alcohol or phenol.
Decanoic acid has the structure CH3(CH2)8C(=O)OH so the acyl fragment found in its ester will be:
CH3(CH2)8C(=O)O (or its deprotonated form).
This classifier looks not only for the decanoate SMARTS fragment but also ensures that the 
methyl group (CH3) is terminal (i.e. not part of a larger alkyl chain).
"""

from rdkit import Chem

def is_decanoate_ester(smiles: str):
    """
    Determines if a molecule is a decanoate ester based on its SMILES string.
    The method first looks for the decanoate acyl fragment via SMARTS and then checks that
    the fragment occurs as a terminal decanoate chain (its CH3 is not extended into a longer chain).
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if a valid decanoate ester substructure is found, False otherwise.
        str: Explanation or reason for the classification.
    """
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS patterns for decanoate ester acyl fragment.
    # The first pattern exactly matches CH3-(CH2)8-C(=O)O
    pattern1 = "[CH3][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2]C(=O)O"
    # The second pattern allows for the ester oxygen being deprotonated (O- instead of O)
    pattern2 = "[CH3][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2]C(=O)[O-]"
    
    smarts_patterns = [pattern1, pattern2]
    
    # Try to find a match for each pattern.
    for pat in smarts_patterns:
        patt = Chem.MolFromSmarts(pat)
        if patt is None:
            continue  # skip invalid pattern
        matches = mol.GetSubstructMatches(patt)
        if not matches:
            continue  # no match for this pattern, try next
        # For each match, check if the matching decanoate moiety is terminal.
        # In our SMARTS, the first atom (index 0) is the CH3 group.
        # For decanoic acid, this CH3 must be terminal (degree = 1 in the acyl group).
        for match in matches:
            ch3_atom = mol.GetAtomWithIdx(match[0])
            # Count neighbors that are not part of this match.
            external_neighbors = [nbr for nbr in ch3_atom.GetNeighbors() if nbr.GetIdx() not in match]
            if len(external_neighbors) == 0:
                # Terminal CH3 group found. This is a valid decanoate moiety.
                return True, "Contains decanoate ester moiety with a terminal decanoate acyl chain"
            # If the CH3 is part of a longer chain, it is not a true decanoate fragment.
    # If no valid match is found, then the decanoate ester fragment is not present properly.
    return False, "Decanoate ester fragment not found or the decanoate acyl group is part of a longer chain"