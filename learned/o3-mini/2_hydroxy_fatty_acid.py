"""
Classifies: CHEBI:10283 2-hydroxy fatty acid
"""
"""
Classifies: Any fatty acid that has a hydroxy functional group in the alpha- or 2-position.
This corresponds to a 2-hydroxy fatty acid.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_2_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is a 2-hydroxy fatty acid based on its SMILES string.
    A 2-hydroxy fatty acid is defined as a fatty acid with a hydroxy group at the alpha carbon
    (i.e. the carbon adjacent to the carboxyl group).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a 2-hydroxy fatty acid, False otherwise.
        str: Reason for the classification decision.
    """
    # Parse the SMILES string to get an RDKit molecule object.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Add explicit hydrogens so that our SMARTS patterns (which depend on explicit H counts) work.
    molH = Chem.AddHs(mol)
    
    # Count the number of carbon atoms (to ensure a minimum chain length typical for fatty acids).
    c_count = sum(1 for atom in molH.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 4:
        return False, "Too few carbons to be considered a fatty acid"
        
    # Define a SMARTS pattern for an alpha-hydroxy carboxylic acid.
    # The pattern looks for a carbon (alpha) that bears an -OH and is directly bonded to a carboxylic acid group.
    alpha_hydroxy_pattern = Chem.MolFromSmarts("[#6]([OX2H])C(=O)[OX1H]")
    if not molH.HasSubstructMatch(alpha_hydroxy_pattern):
        return False, "Missing alpha (2-) hydroxy group adjacent to a carboxyl group"
        
    # Optionally, you can check that the molecule contains a carboxylic acid group.
    # A simple SMARTS for a carboxyl group (using explicit hydrogen on the hydroxyl) is:
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)[OX1H]")
    if not molH.HasSubstructMatch(carboxyl_pattern):
        return False, "Missing a carboxylic acid group"

    # If desired, additional checks such as chain length (e.g., a minimum number of carbons in the alkyl chain)
    # or molecular weight range may be included here to further support the fatty acid classification.

    return True, "Molecule is a 2-hydroxy fatty acid with an alpha-hydroxy group adjacent to the carboxyl group"

# Example: to test the function uncomment the code below
# test_smiles = "CCCCCCCCCCCC[C@@H](O)C(O)=O"  # (R)-2-hydroxytetradecanoic acid
# result, reason = is_2_hydroxy_fatty_acid(test_smiles)
# print(result, reason)