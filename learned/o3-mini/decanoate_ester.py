"""
Classifies: CHEBI:87658 decanoate ester
"""
"""
Classifies: Decanoate ester

A decanoate ester is defined as a fatty acid ester resulting from the condensation 
of decanoic (capric) acid with the hydroxy group of an alcohol or phenol.
Decanoic acid has the structure CH3(CH2)8C(=O)OH so that its ester (decanoate) 
acyl fragment is CH3(CH2)8C(=O)O (or its deprotonated form). This classifier 
searches for that exact substructure with additional checks to ensure that the 
terminal methyl group is truly terminal (using [CH3;D1]) and that the ester oxygen 
indeed bridges to an external group (degree ≥ 2).
"""

from rdkit import Chem

def is_decanoate_ester(smiles: str):
    """
    Determines if a molecule is a decanoate ester based on its SMILES string.
    This is done by looking for a terminal decanoate acyl chain substructure.
    
    The substructure (for both protonated and deprotonated forms) is defined as:
      [CH3;D1][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2]C(=O)O
    or
      [CH3;D1][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2]C(=O)[O-]
      
    The [CH3;D1] restricts the alkyl chain to start with a terminal methyl group.
    Also, after a match is found, we verify that the ester oxygen (the last atom in the match)
    is indeed attached to another group (degree>=2) so that we do not match free decanoic acid.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule contains a decanoate ester moiety as defined, False otherwise.
        str: Explanation or reason for the classification.
    """
    # Parse the SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS patterns for the decanoate acyl fragment.
    # Using [CH3;D1] ensures the methyl is terminal (degree exactly 1).
    pattern_prot = "[CH3;D1][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2]C(=O)O"
    pattern_deprot = "[CH3;D1][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2]C(=O)[O-]"
    
    smarts_patterns = [pattern_prot, pattern_deprot]
    
    for pat in smarts_patterns:
        substruct = Chem.MolFromSmarts(pat)
        if substruct is None:
            continue  # skip invalid pattern (should not happen)
        matches = mol.GetSubstructMatches(substruct)
        if not matches:
            continue  # try next pattern if no match found
        
        # For each match, perform additional checks.
        # Our SMARTS match order is assumed as:
        #  index 0: terminal CH3; indices 1-8: eight CH2's; index 9: carbonyl carbon; index 10: ester oxygen.
        for match in matches:
            if len(match) != 11:
                continue  # unexpected match size
            
            # Check that the terminal CH3 (atom at index match[0]) is indeed terminal.
            ch3_atom = mol.GetAtomWithIdx(match[0])
            if ch3_atom.GetDegree() != 1:
                # It is attached to additional heavy atoms beyond the chain.
                continue
            
            # Check that the ester oxygen (atom at index match[10]) is bridging.
            # It should have at least one neighbor besides the linked acyl carbon.
            oxy_atom = mol.GetAtomWithIdx(match[10])
            if oxy_atom.GetDegree() < 2:
                # This suggests the decanoate chain is not attached to any external R group (e.g. decanoic acid)
                continue

            # If we get here, we've found a valid terminal decanoate ester substructure.
            return True, "Contains decanoate ester moiety with a terminal decanoate acyl chain"
    
    # No valid decanoate ester found.
    return False, "Decanoate ester fragment not found or does not appear as a terminal decanoate acyl chain"