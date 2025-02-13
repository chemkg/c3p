"""
Classifies: CHEBI:35692 dicarboxylic acid
"""
"""
Classifies: Dicarboxylic Acid – Any carboxylic acid containing two carboxy groups.
This improved version uses a SMARTS query to identify a carboxyl group where 
the carbon is double-bonded to an oxygen and singly bonded to an oxygen that is either protonated (–OH) or deprotonated (–O−).
It then applies additional heuristics (e.g. limiting the number of nitrogen atoms, and if no nitrogen is present, requiring a minimum carbon count)
to avoid classifying amino acids or very small acids as dicarboxylic acids.
"""

from rdkit import Chem

def is_dicarboxylic_acid(smiles: str):
    """
    Determines if a molecule qualifies as a dicarboxylic acid (contains exactly two carboxyl groups)
    according to structural heuristics. A carboxyl group is identified as a carbon with a double bond to oxygen
    and a single bond to an oxygen that is either protonated (-OH) or deprotonated (-O−).
    
    Additional heuristics include:
      - If the number of detected carboxyl groups is not exactly 2, the molecule is rejected.
      - If the molecule contains more than 2 nitrogen atoms, it is assumed to be a peptide or complex structure.
      - For molecules without any nitrogen, it is required that there are at least 5 carbon atoms.
      
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule qualifies as a dicarboxylic acid, False otherwise.
        str: An explanation for the classification decision.
    """
    # Parse the SMILES string into a molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern to capture a carboxyl group:
    #   Carbon (C) double bonded to oxygen (=O) and single bonded to an oxygen (O) that carries either a hydrogen or a negative charge.
    carboxyl_smarts = 'C(=O)[O;H1,-]'
    carboxyl_query = Chem.MolFromSmarts(carboxyl_smarts)
    if carboxyl_query is None:
        return False, "Failed to generate SMARTS query for carboxyl group"
    
    # Get all unique substructure matches for the carboxyl group.
    # Each match is a tuple of atom indices corresponding to [C, O] atoms.
    matches = mol.GetSubstructMatches(carboxyl_query, uniquify=True)
    # Use the index of the carbon atom (first in the match tuple) to count unique carboxyl groups.
    carboxyl_carbons = {match[0] for match in matches}
    n_carboxyl = len(carboxyl_carbons)
    
    if n_carboxyl != 2:
        return False, f"Molecule contains {n_carboxyl} carboxyl group(s); need exactly 2."
    
    # Count nitrogen atoms (atomic number 7)
    n_nitrogen = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    # Reject molecules with more than 2 nitrogen atoms (suggesting peptide/complex structure)
    if n_nitrogen > 2:
        return False, "Molecule contains extra nitrogen atoms, suggesting a peptide or complex structure rather than a simple dicarboxylic acid."
    
    # If there are no nitrogen atoms, require at least 5 carbon atoms to avoid very small acids.
    n_carbon = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if n_nitrogen == 0 and n_carbon < 5:
        return False, "Molecule has too few carbon atoms for an interesting dicarboxylic acid."
    
    # All heuristics passed.
    return True, "Molecule contains exactly two carboxyl groups and meets structural heuristics; qualifies as a dicarboxylic acid."

# Example usage (uncomment to test in your own environment):
# print(is_dicarboxylic_acid("OC(=O)C(O)=O"))  # simple dicarboxylic acid test
# print(is_dicarboxylic_acid("C[C@H](NCCC(O)=O)C(O)=O"))  # (S)-beta-alanopine example