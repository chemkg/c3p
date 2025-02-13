"""
Classifies: CHEBI:2571 aliphatic alcohol
"""
"""
Classifies: Aliphatic Alcohol (an alcohol derived from an aliphatic compound)
"""

from rdkit import Chem

def is_aliphatic_alcohol(smiles: str):
    """
    Determines if a molecule is an aliphatic alcohol based on its SMILES string.
    An aliphatic alcohol is defined as an alcohol functional group (-OH) attached
    to an sp3 (non-aromatic) carbon, and not belonging to a carboxyl (acid) group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is classified as an aliphatic alcohol, False otherwise
        str: Explanation for the classification decision
    """
    
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Use SMARTS to find an alcohol group attached to an sp3 (aliphatic) carbon.
    # The pattern "[CX4;!$([CX4](=O))][OX2H]" means:
    #   - [CX4]: a saturated (sp3) carbon atom,
    #   - !$( [CX4](=O) ): the carbon is not part of a carbonyl region (e.g., carboxylic acid),
    #   - [OX2H]: an oxygen atom with two connections and bearing a hydrogen.
    aliphatic_alcohol_pattern = Chem.MolFromSmarts("[CX4;!$([CX4](=O))][OX2H]")
    
    # Check if the molecule contains at least one matching substructure for an aliphatic alcohol.
    if not mol.HasSubstructMatch(aliphatic_alcohol_pattern):
        return False, "No aliphatic alcohol group found (i.e., -OH attached to a non-aromatic sp3 carbon)"
    
    return True, "Molecule contains an aliphatic alcohol group: -OH is attached to an sp3 (non-aromatic) carbon"

# Example usage (you can remove or comment these out when using this as a module):
if __name__ == '__main__':
    # Test with a few example SMILES: methanol (CO) should be classified as aliphatic alcohol
    test_smiles = ["CO", "c1ccccc1O", "CCCCCCCCCCCC(O)CCCC"]  # methanol, phenol, nonadecan-5-ol-like structure
    for smi in test_smiles:
        result, reason = is_aliphatic_alcohol(smi)
        print(f"SMILES: {smi} -> {result} ({reason})")