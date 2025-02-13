"""
Classifies: CHEBI:17984 acyl-CoA
"""
"""
Classifies: Acyl-CoA
Definition: An acyl-CoA is a thioester that results from the formal condensation of 
the thiol group of coenzyme A with the carboxy group of any carboxylic acid.
Example structures include pimeloyl-CoA, stearoyl-CoA, etc.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_acyl_CoA(smiles: str):
    """
    Determines if a molecule is an acyl-CoA based on its SMILES string.
    
    An acyl-CoA is defined as a thioester (R-C(=O)S-) resulting from the condensation
    of the thiol group of coenzyme A with a carboxylic acid.  Its structure should include:
      - A thioester functional group (C(=O)S).
      - A CoA nucleotide moiety—which can be recognized either by the presence of an adenine
        ring (which may appear with various bond orders and substituents) or a part of the pantetheine fragment.
      - Multiple phosphate groups (commonly 2 or more) as part of the nucleotide chain.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as an acyl-CoA, False otherwise.
        str: Explanation for the classification.
    """
    # Parse the SMILES string using RDKit.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 1. Check for the thioester functional group (C(=O)S).
    thioester_smarts = "[CX3](=O)[SX2]"
    thioester = Chem.MolFromSmarts(thioester_smarts)
    if not mol.HasSubstructMatch(thioester):
        return False, "No thioester functional group found"
    
    # 2. Check for the CoA nucleotide part.
    # We'll use several adenine SMARTS patterns to capture different representations.
    adenine_patterns = [
        "n1cnc2c1ncnc2",       # common representation
        "c1ncnc2ncnc12",       # alternative ring closure
        "N1C=NC2=C1N=CN=C2N",   # explicit double bonds
        "c1[nH]c2ncnc(n2)n1"    # sometimes with an explicit [nH] atom
    ]
    has_adenine = any(mol.HasSubstructMatch(Chem.MolFromSmarts(pat)) for pat in adenine_patterns)
    
    # Also try to detect the pantetheine part of CoA.
    # Many acyl-CoA molecules contain the fragment "SCCNC(=O)CCNC(=O)" which is part of pantetheine.
    pantetheine_smarts = "SCCNC(=O)CCNC(=O)"
    pantetheine = Chem.MolFromSmarts(pantetheine_smarts)
    has_pantetheine = mol.HasSubstructMatch(pantetheine)
    
    if not (has_adenine or has_pantetheine):
        return False, "No adenine moiety or pantetheine fragment (CoA nucleotide part) found"
    
    # 3. Check for the presence of phosphate groups.
    # Coenzyme A typically carries several phosphate groups.
    phosphate_smarts = "[OP](=O)(O)"   # matches a phosphate group
    phosphate = Chem.MolFromSmarts(phosphate_smarts)
    phosphate_matches = mol.GetSubstructMatches(phosphate)
    if len(phosphate_matches) < 2:
        return False, "Insufficient phosphate groups for a CoA moiety"
    
    # All key features are present.
    return True, "Molecule contains a thioester linked to a CoA moiety (nucleotide/pantetheine and phosphate groups detected)"

# Example usage (uncomment the lines below for testing):
# test_smiles = "CC(C)(COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12)[C@@H](O)C(=O)NCCC(=O)NCCSC(=O)CCCCCC(O)=O"
# result, reason = is_acyl_CoA(test_smiles)
# print(result, reason)