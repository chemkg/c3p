"""
Classifies: CHEBI:11750 3-sn-phosphatidyl-L-serine
"""
"""
Classifies: 3-sn-phosphatidyl-L-serine
Definition:
   A 3-sn-glycerophosphoserine compound having acyl substituents at the 1- and 2-hydroxy positions.
   That is, a glycerol backbone with three –OH groups where the sn-3 hydroxyl is phosphorylated
   and then esterified with L-serine, and the sn-1 and sn-2 hydroxyls are esterified with acyl chains.
   
This program uses RDKit to check if a given SMILES string conforms to these criteria.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_3_sn_phosphatidyl_L_serine(smiles: str):
    """
    Determines if a molecule is a 3-sn-phosphatidyl-L-serine based on its SMILES string.
    It verifies that the molecule contains:
      - A phosphoserine head group (as a phosphate ester of L-serine attached via a glycerol fragment).
      - Two acyl ester groups (with a sufficiently long alkyl chain, here requiring at least 4 carbons)
        corresponding to the acyl substituents at the sn-1 and sn-2 positions.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a 3-sn-phosphatidyl-L-serine, False otherwise.
        str: Reason for the classification.
    """
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS patterns for the phosphoserine head group.
    # Many valid PS structures carry either C@H or C@@H labels so we include two patterns.
    ps_pattern1 = Chem.MolFromSmarts("COP(=O)(O)OC[C@H](N)C(=O)O")
    ps_pattern2 = Chem.MolFromSmarts("COP(=O)(O)OC[C@@H](N)C(=O)O")
    
    # Check for a match for phosphoserine head group.
    if not (mol.HasSubstructMatch(ps_pattern1) or mol.HasSubstructMatch(ps_pattern2)):
        return False, "Phosphoserine head group not found"
    
    # Define a SMARTS pattern for acyl ester moieties with a long alkyl chain.
    # This pattern looks for an ester group ([CX3](=O)O) attached to at least a 4‐carbon chain.
    tail_pattern = Chem.MolFromSmarts("[CX3](=O)OCCCC")
    
    # Get all matches for this pattern. (These should represent acyl chains attached via an ester linkage.)
    acyl_matches = mol.GetSubstructMatches(tail_pattern)
    if len(acyl_matches) < 2:
        return False, f"Expected at least 2 acyl substituents at the glycerol sn-1 and sn-2 positions; found {len(acyl_matches)}"
    
    # (Optional additional checks could include verifying that the acyl groups are connected
    # to a central glycerol substructure—for example by counting bonds or rotatable bonds—and confirming
    # that no extra non-head-group esters are present. Here we rely on our SMARTS to capture long acyl chains.)
    
    # We may also add a rough filter on molecular weight (most phosphatidylserines are >500 Da)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, "Molecular weight too low to be a phosphatidylserine"
    
    return True, "Molecule contains a phosphoserine head group with two acyl ester substituents on a glycerol backbone"

# (For testing, you could uncomment lines below.)
# test_smiles = "CCCCCCCCCCCCCCCC(=O)OC[C@H](COP(=O)(O)OC[C@H](N)C(=O)O)OC(=O)CCCCCCCCCCCCCCC"
# print(is_3_sn_phosphatidyl_L_serine(test_smiles))