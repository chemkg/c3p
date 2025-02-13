"""
Classifies: CHEBI:17984 acyl-CoA
"""
"""
Classifies: Acyl-CoA
Definition: An acyl-CoA is a thioester that results from the formal condensation 
of the thiol group of coenzyme A with the carboxy group of any carboxylic acid.
Example structures include pimeloyl-CoA, stearoyl-CoA, etc.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_acyl_CoA(smiles: str):
    """
    Determines if a molecule is an acyl-CoA based on its SMILES string.
    
    Acyl-CoA is defined as a thioester (R-C(=O)S-) resulting from the condensation
    of the thiol group of coenzyme A with a carboxylic acid. 
    A valid acyl-CoA should contain:
      - A thioester functional group (C(=O)S).
      - A CoA nucleotide moiety – typically an adenine ring (which may appear with various substituents).
      - Multiple phosphate groups (commonly found in the CoA portion).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as an acyl-CoA, False otherwise.
        str: Explanation for the classification.
    """
    # Parse the SMILES string using RDKit
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 1. Check for the thioester functional group.
    # This pattern finds a carbonyl (C(=O)) directly bonded to a sulfur atom.
    thioester_smarts = "[CX3](=O)[SX2]"
    thioester = Chem.MolFromSmarts(thioester_smarts)
    if not mol.HasSubstructMatch(thioester):
        return False, "No thioester functional group found"
    
    # 2. Check for the CoA nucleotide part.
    # Instead of one strict pattern, we try two adenine (purine) SMARTS to capture various representations.
    adenine_smarts1 = "n1cnc2c(n1)ncnc2"      # a generic pattern for adenine
    adenine_smarts2 = "N1C=NC2=C1N=CN=C2N"    # another representation (often with explicit = bonds)
    adenine1 = Chem.MolFromSmarts(adenine_smarts1)
    adenine2 = Chem.MolFromSmarts(adenine_smarts2)
    
    has_adenine = mol.HasSubstructMatch(adenine1) or mol.HasSubstructMatch(adenine2)
    if not has_adenine:
        return False, "No adenine moiety (CoA nucleotide) found"
    
    # 3. Check for the presence of phosphate groups.
    # Coenzyme A usually carries several phosphate groups.
    phosphate_smarts = "[OP](=O)(O)"   # matches a phosphate group with one double bond and two OH's.
    phosphate = Chem.MolFromSmarts(phosphate_smarts)
    phosphate_matches = mol.GetSubstructMatches(phosphate)
    if len(phosphate_matches) < 2:
        return False, "Insufficient phosphate groups for a CoA moiety"
    
    # If all key features are present, then we classify the molecule as acyl-CoA.
    return True, "Molecule contains a thioester linked to a CoA moiety (adenine nucleotide and phosphate groups detected)"

# Example usage (uncomment the lines below for testing):
# test_smiles = "CC(C)(COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12)[C@@H](O)C(=O)NCCC(=O)NCCSC(=O)CCCCCC(O)=O"
# result, reason = is_acyl_CoA(test_smiles)
# print(result, reason)