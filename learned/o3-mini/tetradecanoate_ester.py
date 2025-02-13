"""
Classifies: CHEBI:87691 tetradecanoate ester
"""
"""
Classifies: tetradecanoate ester

A tetradecanoate ester is defined as a fatty acid ester derived from 
tetradecanoic (myristic) acid. In such a moiety the acid part is a linear, saturated 
chain of exactly 14 carbons (including the carbonyl carbon) that is connected via an ester 
bond to an alcohol (or phenol). This implementation detects an ester substructure and then 
“follows” the acyl chain starting at the carbonyl carbon. Only if that chain is linear (with no 
branching) and exactly 14 carbons long is the molecule classified as containing a tetradecanoate ester.
"""

from rdkit import Chem

def is_tetradecanoate_ester(smiles: str):
    """
    Determines if a molecule (given by its SMILES string) contains an ester 
    group whose acid fragment is exactly tetradecanoic acid (i.e. a linear 14-carbon chain).

    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if a tetradecanoate ester is detected, False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern for an ester group.
    # This pattern captures a carbonyl carbon ([C:1]) attached via a double bond to oxygen and
    # via a single bond to an oxygen ([O:2]). (We will later decide which branch is the acyl chain.)
    ester_smarts = "([C:1](=O)[O:2])"
    ester_query = Chem.MolFromSmarts(ester_smarts)
    if ester_query is None:
        return False, "Error in defining the ester pattern"
    
    # Find all ester matches.
    matches = mol.GetSubstructMatches(ester_query)
    if not matches:
        return False, "No ester functional group found"
    
    # Helper function: from the carbonyl carbon, traverse the connected acyl chain.
    # The traversal starts from a given carbon (the first carbon of the acyl chain) and
    # iteratively moves along a branch provided that:
    #   (1) the bond is single,
    #   (2) the neighbor is a carbon,
    #   (3) the chain does not branch.
    # It returns the number of carbons in this chain (starting count = 1 for the starting atom).
    def count_linear_chain(mol, start_idx, prev_idx):
        count = 1  # count the current atom
        current_idx = start_idx
        coming_from = prev_idx
        while True:
            current_atom = mol.GetAtomWithIdx(current_idx)
            next_neighbors = []
            # Look at bonds on the current atom.
            for bond in current_atom.GetBonds():
                # Only consider single bonds.
                if bond.GetBondType() != Chem.BondType.SINGLE:
                    continue
                other_idx = bond.GetOtherAtomIdx(current_idx)
                # Do not go back where we came from.
                if coming_from is not None and other_idx == coming_from:
                    continue
                other_atom = mol.GetAtomWithIdx(other_idx)
                # Only allow carbon atoms.
                if other_atom.GetAtomicNum() == 6:
                    next_neighbors.append(other_idx)
            # If there is exactly one neighbor, continue along the chain.
            if len(next_neighbors) == 1:
                count += 1
                coming_from = current_idx
                current_idx = next_neighbors[0]
            else:
                # Either no neighbors or a branch is encountered.
                break
        return count

    # Loop over detected ester bonds.
    for match in matches:
        # In our ester SMARTS, match[0] is the carbonyl carbon and match[1] is the ester oxygen.
        carbonyl_idx = match[0]
        ester_oxygen_idx = match[1]
        carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
        
        # Identify the acyl chain branch connected to the carbonyl carbon.
        # The carbonyl typically has three bonds:
        #  - A double bond to a carbonyl oxygen (which is not part of the acyl chain),
        #  - A single bond to the ester oxygen (the other side of the ester),
        #  - And a single bond to the acyl chain carbon.
        acyl_candidates = []
        for nbr in carbonyl_atom.GetNeighbors():
            nbr_idx = nbr.GetIdx()
            # Skip if this neighbor is the ester oxygen.
            if nbr_idx == ester_oxygen_idx:
                continue
            # Also skip if the bond is not single (this will usually skip the double-bonded oxygen).
            bond = mol.GetBondBetweenAtoms(carbonyl_idx, nbr_idx)
            if bond is None or bond.GetBondType() != Chem.BondType.SINGLE:
                continue
            # Only consider carbon atoms.
            if nbr.GetAtomicNum() == 6:
                acyl_candidates.append(nbr_idx)
        
        # We expect exactly one acyl chain carbon attached.
        if len(acyl_candidates) != 1:
            continue
        
        # From the carbonyl carbon, the chain length will be:
        #   1 (for the carbonyl carbon) + length of linear chain starting from the first acyl carbon.
        chain_length = 1 + count_linear_chain(mol, acyl_candidates[0], carbonyl_idx)
        if chain_length == 14:
            return True, ("Contains a tetradecanoate ester moiety "
                          "(detected a linear 14-carbon acyl chain derived from myristic acid)")
    
    # If none of the ester bonds produced a 14-carbon acyl chain, classification fails.
    return False, "No tetradecanoate ester moiety detected"

# For optional testing.
if __name__ == "__main__":
    test_cases = [
        ("C([C@@H](COC(CCCCCCCCCCCCCCC)=O)OC(CCCCCCCCCCCCC)=O)", "DG(16:0/14:0/0:0)"),
        ("CCCCCCCCCCCCCC(=O)OC[C@H](O)COP(O)(=O)OC[C@@H](O)CO", "1,2-ditetradecanoyl-sn-glycerol-3-phosphate"),
        ("CCCCCCCCCCCCCCCC(=O)OCC(CO)OC(=O)CCCCCCCCCCCCCC", "Tetradecanoyl test")
    ]
    for smi, name in test_cases:
        result, reason = is_tetradecanoate_ester(smi)
        print(f"{name}: {result} -> {reason}")