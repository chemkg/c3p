"""
Classifies: CHEBI:86315 methyl sulfide
"""
"""
Classifies: Methyl Sulfide
Definition: Any aliphatic sulfide in which at least one of the organyl groups attached 
to the sulfur is a methyl group.
Note: We exclude molecules that appear to be amino acids or peptides.
"""
from rdkit import Chem

def is_methyl_sulfide(smiles: str):
    """
    Determines if a given molecule is a methyl sulfide.
    
    A methyl sulfide (thioether) by our definition is a molecule that contains at least one
    non‐aromatic sulfur atom (S) that is bonded to two carbon atoms and at least one of them
    is a methyl group (a tetrahedral sp3 carbon bearing exactly three hydrogens). In order to
    reduce false positives from peptides and amino acids (which often have –S–CH3, e.g. L-methionine),
    we first check if the molecule shows an α–amino acid motif.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule qualifies as a methyl sulfide, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens so that hydrogen counts are available.
    mol = Chem.AddHs(mol)
    
    # Exclude molecules that appear to be amino acids or peptides.
    # These patterns look for typical α–amino acid moieties.
    aa_pattern1 = Chem.MolFromSmarts("[C@H](N)C(=O)O")
    aa_pattern2 = Chem.MolFromSmarts("[C@@H](N)C(=O)O")
    if mol.HasSubstructMatch(aa_pattern1) or mol.HasSubstructMatch(aa_pattern2):
        return False, "Molecule appears to be an amino acid or peptide"
    
    # Iterate over sulfur atoms.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 16:
            # Only consider non-aromatic sulfur with exactly two connections.
            if atom.GetIsAromatic():
                continue
            if atom.GetDegree() != 2:
                continue
            
            # For each neighboring atom, check if one is a methyl carbon.
            for neighbor in atom.GetNeighbors():
                # Check neighbor is carbon; use atomic number 6.
                if neighbor.GetAtomicNum() != 6:
                    continue
                # The attached methyl carbon should be tetrahedral (sp3) and non‐aromatic.
                if neighbor.GetHybridization() != Chem.rdchem.HybridizationType.SP3:
                    continue
                if neighbor.GetIsAromatic():
                    continue
                # Count attached hydrogens (from explicit + implicit hydrogens).
                # RDKit has GetTotalNumHs() which sums implicit and explicit hydrogens.
                h_count = neighbor.GetTotalNumHs()
                # Check that this carbon has exactly three hydrogens.
                if h_count == 3:
                    # Found a sulfur atom bonded to a methyl carbon.
                    return True, "Found an aliphatic sulfide with at least one methyl group attached to sulfur."
    
    return False, "No suitable aliphatic sulfide with a methyl group attached was found."
    
# Example usage (uncomment to test):
# test_smiles = "CSCC(=O)[C@H](O)[C@H](O)COP(O)(O)=O"  # 1-(methylthio)ribulose 5-phosphate
# result, reason = is_methyl_sulfide(test_smiles)
# print(result, reason)