"""
Classifies: CHEBI:35366 fatty acid
"""
"""
Classifies: Fatty Acids 
Definition: Any acyclic aliphatic monocarboxylic acid. Natural fatty acids typically have a chain of 4 to 28 carbons (often unbranched and even-numbered), though by extension, all acyclic aliphatic carboxylic acids are included.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_fatty_acid(smiles: str):
    """
    Determines if a molecule is a fatty acid based on its SMILES string.
    
    The molecule must:
      1. Be a valid molecule.
      2. Contain a carboxylic acid functional group (protonated or deprotonated).
      3. Be acyclic (no ring systems).
      4. Possess a carbon chain that is at least 4 carbons long.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule is classified as a fatty acid, False otherwise.
        str: Reason for the classification decision.
    """
    # Attempt to parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern for the carboxylic acid group.
    # This pattern matches a carbon with a double-bonded oxygen and an oxygen that has either one hydrogen or a negative charge.
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)[O;H1,O-]")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No carboxylic acid group found"

    # Check that the molecule is acyclic (no rings).
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule contains ring(s), expected an acyclic fatty acid"
    
    # Count the number of carbon atoms.
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 4:
        return False, f"Too few carbon atoms ({carbon_count}) to be classified as a fatty acid"
    
    # Additional note: while natural fatty acids commonly have 4 to 28 carbons,
    # the definition by extension covers all acyclic aliphatic carboxylic acids.
    # Thus, we do not reject longer chains based on carbon count.
    
    return True, f"Molecule is an acyclic aliphatic carboxylic acid with {carbon_count} carbon(s)"

# Example usage:
if __name__ == "__main__":
    # Test with butyric acid SMILES 
    test_smiles = "CCCC(O)=O"
    result, reason = is_fatty_acid(test_smiles)
    print(result, reason)