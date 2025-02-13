"""
Classifies: CHEBI:10283 2-hydroxy fatty acid
"""
"""
Classifies: Any fatty acid that has a hydroxy functional group in the alpha- or 2-position.
That is, a 2-hydroxy fatty acid.
Examples include 2-hydroxynervonic acid, 2-hydroxytetradecanoic acid, etc.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_2_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is a 2-hydroxy fatty acid based on its SMILES string.
    A 2-hydroxy fatty acid is defined as a fatty acid with a carboxylic acid group
    and an -OH substituent on the carbon alpha (adjacent) to that carboxyl group.
    
    The procedure is:
      1. Parse the SMILES and add explicit hydrogens.
      2. Check basic criteria such as that the molecule is linear and contains enough carbons.
      3. Search for an alpha-hydroxy carboxylic acid substructure.
         Two patterns are attempted:
           a) A tetrahedral (sp3) alpha carbon: [CX4]([OX2H])C(=O)[O;H1]
           b) A trigonal (sp2) alpha carbon: [CX3]([OX2H])C(=O)[O;H1]
      4. If either pattern is found, the molecule is classified as a 2–hydroxy fatty acid.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if the molecule is a 2-hydroxy fatty acid, False otherwise.
        str: Reason for the classification decision.
    """
    # Parse SMILES to create an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens to ensure -OH groups are explicit.
    molH = Chem.AddHs(mol)
    
    # Many fatty acids are linear (acyclic). Reject molecules with rings.
    if molH.GetRingInfo().NumRings() > 0:
        return False, "Molecule contains ring(s), not a typical linear fatty acid"
    
    # Require a minimum number of carbons (adjustable; here we use 4 as a lower bound).
    carbon_count = sum(1 for atom in molH.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 4:
        return False, "Too few carbons to be considered a fatty acid"
    
    # Define SMARTS patterns for the alpha-hydroxy carboxylic acid moiety.
    # Pattern 1: alpha carbon is tetrahedral ([CX4]).
    pattern_sp3 = Chem.MolFromSmarts("[CX4]([OX2H])C(=O)[O;H1]")
    # Pattern 2: alpha carbon is trigonal ([CX3]) as seen in unsaturated examples.
    pattern_sp2 = Chem.MolFromSmarts("[CX3]([OX2H])C(=O)[O;H1]")
    
    # Check if either pattern is present in the molecule.
    if molH.HasSubstructMatch(pattern_sp3) or molH.HasSubstructMatch(pattern_sp2):
        return True, "Molecule contains an alpha-hydroxy group adjacent to a carboxyl group"
    else:
        return False, "Missing alpha (2-) hydroxy group adjacent to a carboxyl group"

# Example usage (uncomment to test):
# test_smiles = "CCCCCCCC\\C=C/CCCCCCCCCCCCC(O)C(O)=O"  # 2-hydroxynervonic acid
# result, reason = is_2_hydroxy_fatty_acid(test_smiles)
# print(result, reason)