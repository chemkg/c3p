"""
Classifies: CHEBI:32957 lysophosphatidic acids
"""
"""
Classifies: lysophosphatidic acids
Definition: Any monoacylglycerol phosphate obtained by hydrolytic removal of one of the two acyl groups of any phosphatidic acid or derivatives therein.
This script uses RDKit to identify a phosphate group, a glycerol backbone fragment,
and exactly one acyl ester (non-phosphate ester) group which defines a monoacyl structure.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_lysophosphatidic_acids(smiles: str):
    """
    Determines if a molecule is a lysophosphatidic acid based on its SMILES string.
    
    The criteria applied are:
    1. The molecule must be valid.
    2. It must contain a phosphate group (e.g., P(=O)(O)(O) moiety).
    3. It must show evidence of a glycerol backbone (i.e. a three-carbon chain bearing hydroxyl groups).
    4. It must have exactly one acyl ester group (an ester not connected to the phosphate group),
       corresponding to the single fatty acyl chain.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule appears to be a lysophosphatidic acid, False otherwise.
        str: Reason for classification.
    """
    
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for the phosphate group.
    # The SMARTS below represents a simple phosphate pattern: P double bonded to O and three single-bonded O's.
    phosphate_smarts = "[P](=O)(O)(O)"
    phosphate_pattern = Chem.MolFromSmarts(phosphate_smarts)
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "Phosphate group not found"

    # Look for a glycerol backbone.
    # Many lysophosphatidic acids contain a glycerol fragment that in simplest form is "OCC(O)CO".
    # This SMARTS pattern is written without chirality so as to match both chiral and non-chiral depictions.
    glycerol_smarts = "OCC(O)CO"
    glycerol_pattern = Chem.MolFromSmarts(glycerol_smarts)
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "Glycerol backbone fragment not found"
    
    # Look for acyl ester group(s) corresponding to the fatty acid chain.
    # We define an acyl ester as an oxygen (which is not attached to a phosphorus) 
    # connected to a carbonyl (C(=O)) that in turn is attached to at least one carbon.
    acyl_smarts = "[O;!$([O]-P)]C(=O)[#6]"
    acyl_pattern = Chem.MolFromSmarts(acyl_smarts)
    acyl_matches = mol.GetSubstructMatches(acyl_pattern)
    n_acyl = len(acyl_matches)
    if n_acyl == 0:
        return False, "No acyl ester group found (expect exactly one fatty acid chain)"
    if n_acyl > 1:
        return False, f"Multiple acyl ester groups found ({n_acyl} groups); expected monoacyl structure"
    
    # Optionally you could add further checks (e.g., sufficient chain length or number of carbons)
    # For instance, many lysophosphatidic acids contain a long chain fatty acid,
    # but we will assume that the one acyl ester group in combination with phosphate and glycerol is enough.
    
    return True, "Molecule has a phosphate group, glycerol backbone, and one acyl ester group characteristic for lysophosphatidic acids"

# Example usage (uncomment to test):
# test_smiles = "P(OC[C@H](O)COC(=O)CCCCCCCCCCCCCCCCCCCCC)(O)(O)=O"  # 1-docosanoyl-glycero-3-phosphate
# result, reason = is_lysophosphatidic_acids(test_smiles)
# print(result, reason)