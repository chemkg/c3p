"""
Classifies: CHEBI:33704 alpha-amino acid
"""
"""
Classifies: Alpha-Amino Acid
Definition: An amino acid in which the amino group is located on the carbon atom at the position alpha to the carboxy group.
This function checks for a carboxyl group (C(=O)[O,OH]) and then looks for a carbon attached to it (alpha carbon) that also has an amino (N) substituent.
"""

from rdkit import Chem

def is_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is an alpha-amino acid based on its SMILES string.
    The criteria used are:
        1. The molecule contains a carboxylic acid moiety (C(=O)[O,OH]).
        2. The carboxyl carbon (of the acid moiety) is attached to a carbon (the alpha carbon),
           which also bears an amino group (N).
    Args:
        smiles (str): SMILES string of the molecule.
    Returns:
        bool: True if the molecule is classified as an alpha-amino acid, False otherwise.
        str: Reason for the classification.
    """
    # Parse SMILES string to RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern for a carboxyl group.
    # The pattern matches a carbon atom double bonded to an oxygen and single bonded to an oxygen (which may be protonated or deprotonated).
    pattern_carboxyl = Chem.MolFromSmarts("C(=O)[O;H1,-]")
    matches = mol.GetSubstructMatches(pattern_carboxyl)
    
    if not matches:
        return False, "No carboxyl group (C(=O)[O,OH]) found in the molecule"
    
    # Check each occurrence of a carboxyl group for the alpha-amino acid pattern.
    for match in matches:
        # In the SMARTS "C(=O)[O;H1,-]", index 0 is the carboxyl carbon.
        carboxyl_idx = match[0]
        carboxyl_atom = mol.GetAtomWithIdx(carboxyl_idx)
        
        # Look for neighbor carbons attached to the carboxyl carbon.
        # The alpha carbon is the one bonded to the carboxyl group.
        alpha_candidates = []
        for neighbor in carboxyl_atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 6:  # must be carbon
                alpha_candidates.append(neighbor)
        
        # For every candidate alpha carbon, check if it carries an amino (nitrogen) group.
        for alpha in alpha_candidates:
            has_amino = False
            for neighbor in alpha.GetNeighbors():
                # Exclude the carboxyl carbon (we already know it)
                if neighbor.GetIdx() == carboxyl_idx:
                    continue
                if neighbor.GetAtomicNum() == 7:  # nitrogen
                    has_amino = True
                    break
            if has_amino:
                return True, ("Alpha amino acid pattern found: "
                              "an alpha carbon (attached to a carboxyl group) also carries an amino group")
    
    return False, "No alpha amino acid pattern found: No alpha carbon with an attached amino group next to a carboxyl group"

# For testing, you can use one of the example SMILES, e.g.:
if __name__ == "__main__":
    test_smiles = [
        "OC(=O)C(N)CN1N=CC=C1",  # 3-(1-Pyrazolyl)-alanine
        "N[C@@H](CC1=CC=C(F)C=C1)C(O)=O",  # 4-fluorophenyl-L-alanine
        "CN[C@@H](Cc1ccccc1)C(O)=O",  # N-methyl-L-phenylalanine
        "OC(=O)CNCC(O)=O"  # iminodiacetic acid, should not match
    ]
    
    for s in test_smiles:
        result, reason = is_alpha_amino_acid(s)
        print(f"SMILES: {s}\nResult: {result}, Reason: {reason}\n")