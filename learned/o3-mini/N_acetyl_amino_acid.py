"""
Classifies: CHEBI:21575 N-acetyl-amino acid
"""
"""
Classifies: N-acetyl-amino acid 
Definition: An N-acyl-amino acid that has acetyl as the acyl group.
The molecule must have an N-acetyl group (CC(=O)N) and an amino acid backbone 
where the acetylated nitrogen is directly attached to an alpha-carbon that bears a carboxylic acid group.
"""

from rdkit import Chem

def is_N_acetyl_amino_acid(smiles: str):
    """
    Determines if a molecule is an N-acetyl-amino acid based on its SMILES string.
    
    The method checks:
    1. The presence of an N-acetyl group (i.e. “CC(=O)N”).
    2. That the nitrogen atom (from the acetyl fragment) is bonded to an alpha carbon
       which, in turn, is bonded to a carboxyl group (-C(=O)O).
       
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is an N-acetyl amino acid, False otherwise.
        str: Explanation of the classification result.
    """
    
    # Parse the SMILES string to a molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # First, look for the N-acetyl moiety.
    # The SMARTS "CC(=O)N" looks for a methyl group (CC) with a carbonyl (C(=O)) attached to a nitrogen.
    acetyl_pattern = Chem.MolFromSmarts("CC(=O)N")
    acetyl_matches = mol.GetSubstructMatches(acetyl_pattern)
    if not acetyl_matches:
        return False, "N-acetyl moiety (CC(=O)N) not found"
    
    # Iterate over each match of the acetyl group.
    # In the match tuple for "CC(=O)N":
    #   index 0: methyl carbon (CH3-)
    #   index 1: carbonyl carbon (C(=O))
    #   index 2: nitrogen (N)
    for match in acetyl_matches:
        n_idx = match[2]
        n_atom = mol.GetAtomWithIdx(n_idx)
        
        # Look for the alpha-carbon attached to the nitrogen.
        # The acetyl nitrogen should be bonded to two atoms: one from the acetyl group (match[1])
        # and one from the amino acid backbone (the alpha-carbon).
        alpha_carbon = None
        for neighbor in n_atom.GetNeighbors():
            # Skip the acetyl carbon.
            if neighbor.GetIdx() == match[1]:
                continue
            # Check that the neighbor is a carbon atom.
            if neighbor.GetAtomicNum() == 6:
                alpha_carbon = neighbor
                break
        if alpha_carbon is None:
            continue  # Try the next occurrence if no alpha-carbon is found.
        
        # Check for a carboxyl group on the alpha-carbon.
        # A carboxylic acid carbon will be attached to two oxygen atoms: one via a double bond (carbonyl)
        # and one via a single bond (hydroxyl or deprotonated oxygen).
        carboxyl_found = False
        o_double_found = False
        o_single_found = False
        
        for neigh in alpha_carbon.GetNeighbors():
            # Look for oxygen atoms.
            if neigh.GetAtomicNum() != 8:
                continue
            bond = mol.GetBondBetweenAtoms(alpha_carbon.GetIdx(), neigh.GetIdx())
            if bond is None:
                continue
            if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                # Potential carbonyl oxygen.
                o_double_found = True
            elif bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                o_single_found = True
        # If both a double-bonded and a single-bonded oxygen are attached, we consider it a carboxyl group.
        if o_double_found and o_single_found:
            carboxyl_found = True
        
        if carboxyl_found:
            return True, "Contains N-acetyl group and amino acid backbone (alpha-carbon with carboxyl group)"
    
    return False, "No alpha-carbon with a carboxyl group linked to an N-acetyl nitrogen was found"
    
# Example usage:
if __name__ == "__main__":
    # Test with one of the provided examples: N-acetyl-L-aspartic acid
    test_smiles = "CC(=O)N[C@@H](CC(O)=O)C(O)=O"
    result, reason = is_N_acetyl_amino_acid(test_smiles)
    print("Result:", result)
    print("Reason:", reason)