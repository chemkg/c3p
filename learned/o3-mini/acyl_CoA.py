"""
Classifies: CHEBI:17984 acyl-CoA
"""
"""
Classifies: Acyl-CoA
Definition: An acyl-CoA is a thioester that results from the formal condensation 
of the thiol group of coenzyme A with the carboxy group of any carboxylic acid.
Example acyl-CoA structures include pimeloyl-CoA, stearoyl-CoA, etc.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_acyl_CoA(smiles: str):
    """
    Determines if a molecule is an acyl-CoA based on its SMILES string.
    
    Acyl-CoA is defined as a thioester (R-C(=O)S-) that results from the condensation of
    the thiol group of coenzyme A with a carboxylic acid. 
    Coenzyme A contains a nucleotide moiety (adenine) and several phosphate groups.
    
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
    
    # Define SMARTS for the thioester group:
    # This pattern finds a carbonyl [CX3](=O) directly bonded to a sulfur atom [SX2]
    thioester_smarts = "[CX3](=O)[SX2]"
    thioester = Chem.MolFromSmarts(thioester_smarts)
    if not mol.HasSubstructMatch(thioester):
        return False, "No thioester functional group found"
    
    # Define SMARTS for an adenine moiety as found in Coenzyme A.
    # Note: many acyl-CoA structures feature the adenine ring (e.g., "n1cnc2c(n1)nc(n2)")
    adenine_smarts = "c1nc2c(n1)nc(n2)"
    adenine = Chem.MolFromSmarts(adenine_smarts)
    if not mol.HasSubstructMatch(adenine):
        return False, "No adenine moiety (CoA nucleotide) found"
        
    # Optionally, check for the presence of phosphate groups common to CoA.
    # Coenzyme A typically contains several phosphate groups; here we require at least two.
    phosphate_smarts = "[OP](=O)(O)"
    phosphate = Chem.MolFromSmarts(phosphate_smarts)
    phos_matches = mol.GetSubstructMatches(phosphate)
    if len(phos_matches) < 2:
        return False, "Insufficient phosphate groups for a CoA moiety"
    
    # If all key features are present, we classify this molecule as an acyl-CoA.
    return True, "Molecule contains a thioester linked to a CoA moiety (adenine and phosphate groups detected)"