"""
Classifies: CHEBI:50760 N-hydroxy-alpha-amino-acid
"""
"""
Classifies: N-hydroxy-alpha-amino-acid
Definition: Any amino acid in which at least one hydrogen attached to the amino group is replaced by a hydroxy group.
Examples include N-hydroxy glycine, N-hydroxy-L-lysine, N,N-dihydroxy-L-isoleucine, etc.
"""
from rdkit import Chem

def is_N_hydroxy_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule (given as a SMILES string) is an N-hydroxy-alpha-amino-acid.
    
    For this classification, the molecule must display an amino acid backbone (an alpha-carbon
    bonded to a carboxylic acid group) and the amino group attached to that alpha carbon must have
    at least one hydroxy (–OH) substituent.
    
    Args:
        smiles (str): SMILES representation of the molecule.
        
    Returns:
        bool: True if the molecule is an N-hydroxy-alpha-amino-acid, False otherwise.
        str: Explanation/reason for the classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."

    # Define SMARTS patterns for the amino acid backbone.
    # We want to capture the alpha carbon and the amino nitrogen via atom-mapping.
    # For most amino acids the alpha carbon has one hydrogen (chiral center):
    aa_pattern_non_gly = Chem.MolFromSmarts("[C:1;H1]([#7:2])C(=O)[O]")
    # Glycine (achiral) has two hydrogens on the alpha carbon:
    aa_pattern_gly = Chem.MolFromSmarts("[C:1;H2]([#7:2])C(=O)[O]")
    
    # Try to find a match with either pattern
    matches = mol.GetSubstructMatches(aa_pattern_non_gly)
    if not matches:
        matches = mol.GetSubstructMatches(aa_pattern_gly)
    if not matches:
        return False, "No amino acid backbone detected (missing alpha-carbon with adjacent carboxyl group and amino group)."
    
    # Loop through each match; the second atom in the match is the amino nitrogen.
    for match in matches:
        alpha_carbon_idx = match[0]
        amino_nitrogen_idx = match[1]
        amino_nitrogen = mol.GetAtomWithIdx(amino_nitrogen_idx)
        
        # Check if the amino nitrogen has at least one oxygen substituent (via a single bond).
        # This is our proxy for an -OH substitution on the amino group.
        for neighbor in amino_nitrogen.GetNeighbors():
            # Skip if the neighbor is the alpha carbon (or any carbon in the backbone)
            if neighbor.GetIdx() == alpha_carbon_idx:
                continue
            if neighbor.GetAtomicNum() == 8:  # oxygen
                # Get the bond between the nitrogen and the oxygen
                bond = mol.GetBondBetweenAtoms(amino_nitrogen_idx, neighbor.GetIdx())
                if bond is not None and bond.GetBondType() == Chem.BondType.SINGLE:
                    # Check if the oxygen carries at least one hydrogen.
                    # (RDKit uses implicit Hs so GetTotalNumHs() works for our purposes.)
                    if neighbor.GetTotalNumHs() > 0:
                        return True, ("Molecule contains an amino acid backbone and the amino " 
                                      "group has a hydroxy substituent (N-hydroxy modification detected).")
    return False, "Amino acid backbone detected but no N-hydroxy substitution on the amino group found."

# Example usage:
if __name__ == "__main__":
    test_smiles = [
        "ONCC(O)=O",  # N-hydroxyglycine
        "N[C@@H](CCCCNO)C(O)=O",  # N(6)-hydroxy-L-lysine
        "CC(C)[C@H](N(O)O)C(O)=O",  # N,N-dihydroxy-L-valine
        "N[C@@H](CCCC)C(O)=O"  # Regular amino acid lacking N-hydroxy substituent
    ]
    for smi in test_smiles:
        result, reason = is_N_hydroxy_alpha_amino_acid(smi)
        print(f"SMILES: {smi}\nResult: {result}\nReason: {reason}\n")