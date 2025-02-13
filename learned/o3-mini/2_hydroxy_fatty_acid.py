"""
Classifies: CHEBI:10283 2-hydroxy fatty acid
"""
"""
Classifies: Any fatty acid that has a hydroxy functional group in the alpha- or 2-position,
i.e. a 2-hydroxy fatty acid.
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

    # Add explicit hydrogens so that SMARTS patterns work with explicit hydrogen counts.
    molH = Chem.AddHs(mol)
    
    # Make sure there is a minimum chain length typical for fatty acids.
    c_count = sum(1 for atom in molH.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 4:
        return False, "Too few carbons to be considered a fatty acid"

    # Define a SMARTS pattern for the alpha-hydroxy fatty acid motif.
    # This pattern looks for a CH group (alpha carbon with one hydrogen) that bears an -OH
    # and is directly bonded to a carboxyl group with a carbonyl and an oxygen that is either protonated or negative.
    alpha_hydroxy_pattern = Chem.MolFromSmarts("[CH]([OX2H])C(=O)[O;H1,O-]")
    if not molH.HasSubstructMatch(alpha_hydroxy_pattern):
        return False, "Missing alpha (2-) hydroxy group adjacent to a carboxyl group"
    
    # Additionally, check that a carboxylic acid group is present.
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)[O;H1,O-]")
    if not molH.HasSubstructMatch(carboxyl_pattern):
        return False, "Missing a carboxylic acid group"

    # Optional: Check that the molecular weight and alkyl chain are in the fatty acid range.
    mol_wt = rdMolDescriptors.CalcExactMolWt(molH)
    if mol_wt < 100:  # most fatty acids are heavier than about 100 Da
        return False, "Molecular weight too low for a fatty acid"
    
    return True, "Molecule is a 2-hydroxy fatty acid with an alpha-hydroxy group adjacent to the carboxyl group"

# Example usage (uncomment the lines below to test):
# test_smiles = "CCCCCCCC\\C=C/CCCCCCCCCCCCC(O)C(O)=O"  # 2-hydroxynervonic acid
# result, reason = is_2_hydroxy_fatty_acid(test_smiles)
# print(result, reason)