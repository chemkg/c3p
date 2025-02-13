"""
Classifies: CHEBI:48030 tetrapeptide
"""
"""
Classifies: Tetrapeptide
Definition: Any molecule that contains four amino‐acid residues connected by peptide linkages.
A linear tetrapeptide should have three peptide bonds in its main chain.
"""

from rdkit import Chem

def is_tetrapeptide(smiles: str):
    """
    Determines if the given molecule (as SMILES) is a tetrapeptide.
    The strategy is to look for backbone peptide bonds. In most linear tetrapeptides,
    there are three such bonds (connecting four amino acid residues).
    
    Args:
        smiles (str): The SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is a tetrapeptide, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a SMARTS pattern for a backbone peptide linkage.
    # The pattern "N-C(=O)-[C]" searches for an amide bond (peptide bond)
    # where the nitrogen is connected to a carbonyl carbon and then to a (typically alpha) carbon.
    peptide_linkage = Chem.MolFromSmarts("N-C(=O)-[C]")
    if peptide_linkage is None:
        return False, "Could not construct peptide linkage SMARTS pattern"
    
    # Get all substructure matches for the backbone pattern
    matches = mol.GetSubstructMatches(peptide_linkage)
    if not matches:
        return False, "No peptide bonds found"
    
    # Each match is a triple (N_idx, C_idx, Ca_idx) corresponding to N-C(=O)-C.
    # In a standard tetrapeptide, there should be three peptide bonds connecting four residues.
    # (Note: if the molecule contains additional amide bonds from side chain groups, this heuristic could be confounded.)
    peptide_bonds = set()
    for match in matches:
        # Use the pair (carbonyl carbon, peptide nitrogen) to represent the peptide bond.
        # The peptide bond is between the carbonyl carbon (match[1]) and the nitrogen (match[0]).
        peptide_bonds.add((match[1], match[0]))
    
    n_peptide_bonds = len(peptide_bonds)
    if n_peptide_bonds != 3:
        return False, f"Found {n_peptide_bonds} backbone peptide bond(s); expected exactly 3 for a tetrapeptide"
    
    # Optionally, one might also check that the molecule’s overall size is within the expected range.
    # For many tetrapeptides the molecular weight is within a moderate range.
    # Here we skip that check.
    
    return True, "Molecule contains four amino-acid residues connected by 3 peptide bonds (backbone) consistent with a tetrapeptide"

# Example usage:
if __name__ == "__main__":
    # Example tetrapeptide (Glu-Lys-Trp-Ala)
    example_smiles = "C[C@H](NC(=O)[C@H](Cc1c[nH]c2ccccc12)NC(=O)[C@H](CCCCN)NC(=O)[C@@H](N)CCC(O)=O)C(O)=O"
    result, reason = is_tetrapeptide(example_smiles)
    print(result, reason)