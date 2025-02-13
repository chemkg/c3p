"""
Classifies: CHEBI:87691 tetradecanoate ester
"""
"""
Classifies: tetradecanoate ester
A tetradecanoate ester is defined as a fatty acid ester obtained by condensation
of the carboxy group of tetradecanoic acid (myristic acid) with a hydroxy group
of an alcohol or phenol.

This implementation first locates ester substructures (using a SMARTS pattern to
identify a carbonyl carbon linked to an ester oxygen) and then, for each candidate,
it “follows” the acyl chain starting at the carbonyl carbon. The algorithm accepts
the match only if the acyl chain is a linear, saturated alkyl chain that has exactly
14 carbons (the carbonyl carbon plus 13 other sp3 carbons, representing CH3-(CH2)12-C(=O)O).
This helps avoid false positives where a longer chain (or a fragment of a longer chain)
would otherwise match.
"""

from rdkit import Chem

def is_tetradecanoate_ester(smiles: str):
    """
    Determines if a molecule contains an ester substructure where the acid part 
    is exactly tetradecanoic acid (myristic acid), i.e. it comprises a linear, saturated
    14-carbon chain (including the carbonyl carbon) attached via an ester oxygen.

    This function locates ester bonds defined by a SMARTS pattern and then, for each
    candidate ester, it tries to traverse (linearly) the acyl chain. If exactly 14 carbons
    (the carbonyl and 13 additional carbons) are present in that chain, we classify the
    molecule as containing a tetradecanoate ester.

    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if a tetradecanoate ester moiety is detected, False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # SMARTS pattern for an ester group: the pattern captures a carbonyl carbon ([C:1])
    # that is double-bonded to an oxygen and single-bonded to an oxygen ([O:2]).
    ester_smarts = "([C:1](=O)[O:2])"
    ester_query = Chem.MolFromSmarts(ester_smarts)
    if ester_query is None:
        return False, "Error in defining the ester pattern"
    
    matches = mol.GetSubstructMatches(ester_query)
    if not matches:
        return False, "No ester functional group found"
    
    # Helper function: recursively traverse a linear chain of sp3 carbon atoms.
    # It returns the number of carbon atoms in this chain including the current atom.
    def get_linear_chain_length(atom, prev):
        """
        Traverse the acyl chain starting from a given carbon atom.
        Only if the connected chain is linear (no branching) and consists solely
        of saturated carbon atoms are counts returned.
        
        Args:
            atom: current RDKit atom.
            prev: the previous atom (to avoid going backwards).
            
        Returns:
            int: number of carbon atoms in the linear chain (including current) if linear,
                 otherwise None if branching or non-carbon is encountered.
        """
        # Only count carbon atoms.
        if atom.GetAtomicNum() != 6:
            return None
        
        # Get neighboring carbons (exclude the atom we came from).
        neighbors = [nbr for nbr in atom.GetNeighbors() 
                     if nbr.GetAtomicNum() == 6 and (prev is None or nbr.GetIdx() != prev.GetIdx())]
        # Allow terminal: if no further carbon neighbor, count this atom.
        if len(neighbors) == 0:
            return 1
        # If more than one carbon neighbor, the chain is branched.
        if len(neighbors) > 1:
            return None
        # For the single next carbon, ensure the bond is a single bond.
        bond = mol.GetBondBetweenAtoms(atom.GetIdx(), neighbors[0].GetIdx())
        if bond is None or bond.GetBondType() != Chem.BondType.SINGLE:
            return None
        # Recursively count.
        next_count = get_linear_chain_length(neighbors[0], atom)
        if next_count is None:
            return None
        return 1 + next_count

    # Loop over each ester match to see if any gives a linear chain of exactly 14 carbons.
    # In a proper tetradecanoate, the acid fragment is CH3-(CH2)12-C(=O) (14 carbons total).
    for match in matches:
        # match[0] corresponds to the carbonyl carbon and match[1] to the ester oxygen.
        carbonyl_atom = mol.GetAtomWithIdx(match[0])
        ester_oxygen = mol.GetAtomWithIdx(match[1])
        
        # From the carbonyl atom, identify the acyl chain neighbor.
        # The carbonyl carbon is typically bonded to:
        #   - a double-bonded oxygen (not part of the acyl chain)
        #   - the ester oxygen (already used)
        #   - the acyl chain carbon.
        chain_neighbors = []
        for nbr in carbonyl_atom.GetNeighbors():
            # Skip if neighbor is the ester oxygen (or if the bond is double).
            bond = mol.GetBondBetweenAtoms(carbonyl_atom.GetIdx(), nbr.GetIdx())
            if nbr.GetIdx() == ester_oxygen.GetIdx():
                continue
            # Also skip if the bond is a double bond (the carbonyl double-bonded oxygen).
            if bond and bond.GetBondType() != Chem.BondType.SINGLE:
                continue
            if nbr.GetAtomicNum() == 6:  # we want a carbon.
                chain_neighbors.append(nbr)
        
        # We expect exactly one acyl chain connected.
        if len(chain_neighbors) != 1:
            continue
        
        # Perform a depth-first, linear traversal starting at the carbonyl atom.
        # We want to count the number of carbon atoms in the acyl fragment,
        # including the carbonyl carbon.
        chain_length = get_linear_chain_length(carbonyl_atom, None)
        if chain_length == 14:
            return True, "Contains a tetradecanoate ester moiety (myristate ester group detected)"
    
    return False, "No tetradecanoate ester moiety detected"

# Optional: allow some testing when run as a script.
if __name__ == "__main__":
    # Example test cases:
    test_cases = [
        ("C([C@@H](COC(CCCCCCCCCCCCCCC)=O)OC(CCCCCCCCCCCCC)=O)O", "DG(16:0/14:0/0:0)"),
        ("CCCCCCCCCCCCCC(=O)OC[C@H](O)COP(O)(=O)OC[C@@H](O)CO", "1-tetradecanoyl-sn-glycero-3-phospho-(1'-sn-glycerol)"),
        ("O(C(=O)CCCCCCCCCCCCCCC)C[C@@H](OC(=O)CCCCCCC/C=C\\CCCCCCCC)COC(=O)CCCCCCC/C=C\\C/C=C\\CCCCC", "TG(16:0/18:1(9Z)/18:2(9Z,12Z))[iso6] (false positive example)")
    ]
    for smi, name in test_cases:
        result, reason = is_tetradecanoate_ester(smi)
        print(f"{name}: {result} -> {reason}")