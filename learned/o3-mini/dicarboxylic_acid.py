"""
Classifies: CHEBI:35692 dicarboxylic acid
"""
"""
Classifies: Dicarboxylic Acid – Any carboxylic acid containing two carboxy groups.
Improved approach:
  • Uses an expanded SMARTS that matches both protonated and deprotonated carboxyl groups.
  • Applies an extra heuristic filter on the number of nitrogen atoms (peptide-like molecules
    often have more than 2 nitrogens) in order to reduce false positives.
"""
from rdkit import Chem

def is_dicarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a dicarboxylic acid based on its SMILES string.
    A dicarboxylic acid is defined for our purposes as any molecule containing exactly two
    carboxyl groups (–COOH or –COO–) and not being part of a larger peptide-like structure.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is classified as a dicarboxylic acid, False otherwise.
        str: Explanation of the classification decision.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern for a carboxyl group.
    # This pattern matches [CX3](=O) followed by an oxygen that is either protonated (-OH) or deprotonated (-O-).
    carboxyl_smarts = "[CX3](=O)[OX1H0-,OX2H1]"
    carboxyl_group = Chem.MolFromSmarts(carboxyl_smarts)
    if carboxyl_group is None:
        return False, "Error creating SMARTS for carboxyl group"
    
    # Get the matches for the carboxyl group.
    matches = mol.GetSubstructMatches(carboxyl_group)
    n_carboxyl = len(matches)
    
    # Heuristic: count nitrogen atoms. Many peptides and biomolecules tend to have >2 nitrogens.
    n_nitrogen = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    
    # Classification:
    # First, require that exactly 2 carboxyl groups are present.
    if n_carboxyl != 2:
        return False, f"Molecule contains {n_carboxyl} carboxyl group(s); need exactly 2."
    
    # Next, exclude molecule that are likely peptide-like.
    if n_nitrogen > 2:
        return False, "Molecule contains extra nitrogen atoms, suggesting a peptide or complex structure rather than a simple dicarboxylic acid."
    
    # If both conditions are met, it qualifies.
    return True, "Molecule contains exactly two carboxyl groups and has limited nitrogen content; qualifies as a dicarboxylic acid."

# Example usage:
# Uncomment the lines below to try out a couple of examples in your environment.
# print(is_dicarboxylic_acid("OC(=O)C(O)=O"))  # e.g., simple dicarboxylic acid
# print(is_dicarboxylic_acid("C[C@H](NCCC(O)=O)C(O)=O"))  # (S)-beta-alanopine, should return True