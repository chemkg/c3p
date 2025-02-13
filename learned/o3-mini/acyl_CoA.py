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
    of the thiol group of coenzyme A with a carboxylic acid. Its structure should include:
      - A thioester functional group (C(=O)S).
      - A CoA nucleotide part, which may be recognized either by the presence of an adenine
        ring (in varying representations) or by a pantetheine fragment.
      - Multiple phosphate groups (typically at least three phosphorus atoms in the CoA moiety).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as an acyl-CoA, False otherwise.
        str: Explanation for the classification.
    """
    # Parse the SMILES using RDKit.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 1. Check for the thioester functional group: a carbonyl directly bound to a sulfur.
    thioester_smarts = "[CX3](=O)[SX2]"
    thioester = Chem.MolFromSmarts(thioester_smarts)
    if not mol.HasSubstructMatch(thioester):
        return False, "No thioester functional group found"
    
    # 2. Check for the CoA nucleotide part.
    # Try to capture the adenine ring with several SMARTS patterns.
    adenine_patterns = [
        "n1cnc2c1ncnc2",      # common adenine representation
        "c1ncnc2ncnc12", 
        "N1C=NC2=C1N=CN=C2N",  
        "c1[nH]c2ncnc(n2)n1"   
    ]
    has_adenine = any(mol.HasSubstructMatch(Chem.MolFromSmarts(pat)) for pat in adenine_patterns)
    
    # Recognize part of the pantetheine fragment which is also common in CoA.
    pantetheine_smarts = "SCCNC(=O)CCNC(=O)"
    pantetheine = Chem.MolFromSmarts(pantetheine_smarts)
    has_pantetheine = mol.HasSubstructMatch(pantetheine)
    
    if not (has_adenine or has_pantetheine):
        return False, "No adenine moiety or pantetheine fragment (CoA nucleotide part) found"
    
    # 3. Check for the presence of phosphate groups.
    # Rather than matching a rigid phosphate SMARTS pattern, we count phosphorus atoms.
    phosphate_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
    if phosphate_count < 3:
        return False, "Insufficient phosphate groups for a CoA moiety (found {} phosphorus atom(s))".format(phosphate_count)
    
    # If all three key features are detected then classify as acyl-CoA.
    return True, "Molecule contains a thioester linked to a CoA moiety (CoA nucleotide part and sufficient phosphate groups detected)"

# Example usage:
# Uncomment the following lines to test with a known acyl-CoA molecule (pimeloyl-CoA example)
# test_smiles = "CC(C)(COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12)[C@@H](O)C(=O)NCCC(=O)NCCSC(=O)CCCCCC(O)=O"
# result, reason = is_acyl_CoA(test_smiles)
# print(result, reason)