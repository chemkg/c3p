"""
Classifies: CHEBI:17962 CDP-diacylglycerol
"""
"""
Classifies: CDP-diacylglycerol
Definition: A CDP-glycerol having unspecified acyl groups (most commonly fatty acyl groups) at the 1- and 2-positions.
Heuristic rules used:
  1. The molecule must be parsable.
  2. It must contain a cytidine moiety (using the SMARTS "n1ccc(N)nc1=O").
  3. It must have a bridging diphosphate (P–O–P) motif (SMARTS "[P]-O-[P]").
  4. It must contain at least two ester bonds of the type acyl group (SMARTS "OC(=O)[#6]") that are assumed to represent the fatty acyl groups.
"""

from rdkit import Chem

def is_CDP_diacylglycerol(smiles: str):
    """
    Determines if a molecule is a CDP-diacylglycerol based on its SMILES string.
    
    A CDP-diacylglycerol must contain a cytidine diphosphate head group (i.e. a cytidine
    moiety and a diphosphate bridge) and a glycerol backbone substituted at the 1- and 2-positions
    with fatty acyl (ester) groups.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a CDP-diacylglycerol, otherwise False.
        str: Reason for the classification decision.
    """
    # Parse SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for the presence of a cytidine moiety.
    # This SMARTS picks up the cytosine base component as found in many CDP head groups.
    cytidine_pattern = Chem.MolFromSmarts("n1ccc(N)nc1=O")
    if not mol.HasSubstructMatch(cytidine_pattern):
        return False, "Cytidine moiety not found"
    
    # Check for the diphosphate bridge by searching for a P–O–P connectivity.
    pp_bridge_pattern = Chem.MolFromSmarts("[P]-O-[P]")
    if not mol.HasSubstructMatch(pp_bridge_pattern):
        return False, "Diphosphate (P–O–P) bridge not found"
    
    # Count acyl (ester) groups that are attached to glycerol.
    # We look for the ester moiety pattern that is common to fatty acyl chains: O-C(=O)-[C]
    acyl_pattern = Chem.MolFromSmarts("OC(=O)[#6]")
    acyl_matches = mol.GetSubstructMatches(acyl_pattern)
    if len(acyl_matches) < 2:
        return False, f"Found {len(acyl_matches)} acyl ester group(s); need at least 2 for diacylglycerol"
    
    # If all checks pass, we conclude it is a CDP-diacylglycerol.
    return True, "Contains cytidine headgroup with diphosphate bridge and at least 2 acyl ester groups"
    
# The function above can be tested with one of the provided example SMILES.
if __name__ == "__main__":
    # Example: CDP-DG(16:0/18:2(9Z,12Z))
    smiles_test = "P(OC[C@H]1O[C@@H](N2C=CC(=NC2=O)N)C(O)[C@H]1O)(OP(OC[C@H](OC(=O)CCCCCCC/C=C\\C/C=C\\CCCCC)COC(=O)CCCCCCCCCCCCCCC)(O)=O)(O)=O"
    result, reason = is_CDP_diacylglycerol(smiles_test)
    print(result, reason)