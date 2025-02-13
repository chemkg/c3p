"""
Classifies: CHEBI:35366 fatty acid
"""
"""
Classifies: Fatty Acids 
Definition: Any acyclic aliphatic carboxylic acid. Although natural fatty acids are typically monocarboxylic,
they may appear in esterified form and a few derivatives (e.g. hydroxylated, aminated) are tolerated.
This classifier ensures that the molecule is:
  • A valid molecule.
  • Acyclic.
  • Contains exactly one terminal carboxylic acid group.
  • Does not contain amide bonds (to avoid peptides).
  • Has a sufficient number of carbon atoms and a carbon-to-heteroatom ratio indicative of a fatty acid.
Examples of accepted fatty acids (with various chain lengths) are included in the outcomes.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_fatty_acid(smiles: str):
    """
    Determines if a molecule is a fatty acid based on its SMILES string.
    
    The molecule must:
      1. Be a valid molecule.
      2. Be acyclic (no rings).
      3. Contain a terminal carboxylic acid group – a carbon with a double-bonded oxygen and an -OH (or O-)
         which is attached to exactly one carbon (i.e. at the end of a chain).
      4. Not contain amide bonds (e.g. C(=O)N) which are indicative of peptide bonds.
      5. Contain a sufficiently long carbon chain (at least 4 carbons) and have a carbon-to-heteroatom ratio 
         that is typical for an aliphatic acid.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a fatty acid, False otherwise.
        str: Reason for the classification decision.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check that the molecule is acyclic.
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule contains ring(s), expected an acyclic fatty acid"
    
    # Define SMARTS for a carboxylic acid group.
    # This pattern matches a carbon with a double-bonded oxygen and a singly-bonded oxygen (with hydrogen or negative charge).
    carboxyl_smarts = "C(=O)[O;H1,O-]"
    carboxyl_query = Chem.MolFromSmarts(carboxyl_smarts)
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_query)
    
    # Ensure that the carboxyl group is terminal.
    # For each matching carboxyl carbon (match[0] is the carbon), count the number of carbon neighbors.
    terminal_carboxyl_matches = []
    for match in carboxyl_matches:
        c_idx = match[0]
        c_atom = mol.GetAtomWithIdx(c_idx)
        # Count how many neighbors of the carboxyl carbon are carbon atoms.
        carbon_neighbors = [nbr for nbr in c_atom.GetNeighbors() if nbr.GetAtomicNum() == 6]
        if len(carbon_neighbors) == 1:
            terminal_carboxyl_matches.append(match)
    
    if len(terminal_carboxyl_matches) != 1:
        return False, f"Found {len(terminal_carboxyl_matches)} terminal carboxyl group(s), expected exactly 1"
    
    # Reject molecules with amide bonds, as these are common in peptides and proteins.
    amide_query = Chem.MolFromSmarts("C(=O)N")
    if mol.HasSubstructMatch(amide_query):
        return False, "Molecule contains amide bond(s), likely not a fatty acid"
    
    # Count the number of carbon atoms.
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 4:
        return False, f"Too few carbon atoms ({carbon_count}) to be a fatty acid"
    
    # As a heuristic for “aliphatic character”, compute the ratio of carbon atoms to heteroatoms 
    # (exclude hydrogen). Fatty acids usually have many more carbons than heteroatoms.
    hetero_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() not in (1, 6))
    ratio = carbon_count / hetero_count if hetero_count > 0 else float('inf')
    if ratio < 1.5:
        return False, f"Low carbon to heteroatom ratio ({ratio:.2f}), not typical for a fatty acid"
    
    return True, f"Molecule is an acyclic aliphatic carboxylic acid with {carbon_count} carbon(s)"

# Example usage:
if __name__ == "__main__":
    # Test with butyric acid SMILES
    test_smiles = "CCCC(O)=O"
    result, reason = is_fatty_acid(test_smiles)
    print(result, reason)