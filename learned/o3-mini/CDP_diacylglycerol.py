"""
Classifies: CHEBI:17962 CDP-diacylglycerol
"""
"""
Classifies: CDP-diacylglycerol
Definition: A CDP-glycerol having unspecified acyl groups (most commonly fatty acyl groups) 
at the 1- and 2-positions.
Heuristic improvements in this implementation:
  1. The molecule must be parsable.
  2. It must contain a cytidine moiety (SMARTS "n1ccc(N)nc1=O").
  3. It must have a diphosphate bridge. Instead of using the loose pattern "[P]-O-[P]",
     we use a more specific pattern "P(=O)(O)[O]-[P](=O)(O)" which ensures correct phosphoryl groups.
  4. It must show evidence for a glycerol attachment from the headgroup. Rather than "OP(OC)",
     we now require a branch containing an "OCC" fragment (as in "-OP(OCC)"), a rough indicator for glycerol.
  5. It must contain at least two acyl ester groups (SMARTS "OC(=O)[#6]"), representing the fatty acyl groups.
"""

from rdkit import Chem

def is_CDP_diacylglycerol(smiles: str):
    """
    Determines if a molecule is a CDP-diacylglycerol based on its SMILES string.
    
    A CDP-diacylglycerol must contain:
      - A cytidine moiety,
      - A diphosphate bridge with the correct phosphoryl pattern,
      - A glycerol attachment (indicated by a branch with an OCC fragment),
      - And at least two acyl ester groups.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a CDP-diacylglycerol, otherwise False.
        str: Reason for the classification decision.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for the cytidine moiety.
    # This pattern picks up the cytosine base component that is present in the CDP headgroup.
    cytidine_pattern = Chem.MolFromSmarts("n1ccc(N)nc1=O")
    if not mol.HasSubstructMatch(cytidine_pattern):
        return False, "Cytidine moiety not found"
    
    # Check for the diphosphate bridge.
    # We now use a more chemically specific pattern that requires phosphoryl groups:
    # pattern: P(=O)(O)[O]-[P](=O)(O)
    pp_bridge_pattern = Chem.MolFromSmarts("P(=O)(O)[O]-[P](=O)(O)")
    pp_matches = mol.GetSubstructMatches(pp_bridge_pattern)
    if not pp_matches:
        return False, "Diphosphate (P(=O)(O)[O]-[P](=O)(O)) bridge not found"
    if len(pp_matches) != 1:
        return False, f"Expected exactly one diphosphate bridge but found {len(pp_matches)}"
    
    # Check for a glycerol attachment indicator.
    # Instead of a loose "OP(OC)", we require that one branch of the phosphate shows "OCC"
    # which is more indicative that the headgroup is linked to a glycerol (CH2CH(OH)…) moiety.
    glycerol_pattern = Chem.MolFromSmarts("OP(OCC)")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "Glycerol attachment not detected (missing typical OP(OCC) fragment)"
    
    # Count the acyl ester groups.
    # We look for groups of the form OC(=O)[#6] which represent the fatty acyl esters.
    acyl_pattern = Chem.MolFromSmarts("OC(=O)[#6]")
    acyl_matches = mol.GetSubstructMatches(acyl_pattern)
    if len(acyl_matches) < 2:
        return False, f"Found {len(acyl_matches)} acyl ester group(s); need at least 2"
    
    return True, "Contains cytidine headgroup with diphosphate bridge, glycerol attachment, and at least 2 acyl ester groups"

# Testing with one example (CDP-DG(16:0/18:2(9Z,12Z))):
if __name__ == "__main__":
    test_smiles = "P(OC[C@H]1O[C@@H](N2C=CC(=NC2=O)N)C(O)[C@H]1O)(OP(OC[C@H](OC(=O)CCCCCCC/C=C\\C/C=C\\CCCCC)COC(=O)CCCCCCCCCCCCCCC)(O)=O)(O)=O"
    result, reason = is_CDP_diacylglycerol(test_smiles)
    print(result, reason)