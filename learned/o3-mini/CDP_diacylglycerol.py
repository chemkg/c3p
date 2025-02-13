"""
Classifies: CHEBI:17962 CDP-diacylglycerol
"""
"""
Classifies: CDP-diacylglycerol
Definition: A CDP-glycerol having unspecified acyl groups (most commonly fatty acyl groups)
at the 1- and 2-positions.
This classifier ensures that the molecule contains:
  - A cytidine moiety,
  - A diphosphate bridge (with exactly one unique group of two phosphorus atoms connected via an oxygen),
  - A glycerol attachment indicator (a branch with an OCC fragment),
  - And at least two acyl ester groups.
"""

from rdkit import Chem

def is_CDP_diacylglycerol(smiles: str):
    """
    Determines if a molecule is a CDP-diacylglycerol based on its SMILES string.
    
    A CDP-diacylglycerol must contain:
      - A cytidine moiety,
      - A diphosphate bridge, checked using a specific SMARTS pattern and deduplication
        (exactly one unique set of two phosphorus atoms connected via an oxygen),
      - A glycerol attachment indicator (e.g., an OP(OCC) fragment),
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
    
    # Check for the cytidine moiety in the headgroup.
    cytidine_pattern = Chem.MolFromSmarts("n1ccc(N)nc1=O")
    if not mol.HasSubstructMatch(cytidine_pattern):
        return False, "Cytidine moiety not found"
    
    # Identify the diphosphate bridge.
    # Use a SMARTS that specifies two phosphorus atoms with phosphoryl groups joined by an oxygen.
    pp_bridge_pattern = Chem.MolFromSmarts("P(=O)(O)[O]-[P](=O)(O)")
    pp_matches = mol.GetSubstructMatches(pp_bridge_pattern)
    if not pp_matches:
        return False, "Diphosphate (P(=O)(O)[O]-[P](=O)(O)) bridge not found"
    
    # Deduplicate to find unique diphosphate bridges.
    unique_pp_groups = set()
    for match in pp_matches:
        # The match returns a tuple of two atom indices (for the two phosphorus atoms in the pattern).
        # We create a frozenset so that (a,b) and (b,a) are considered the same.
        unique_pp_groups.add(frozenset(match))
    
    if len(unique_pp_groups) != 1:
        return False, f"Expected exactly one diphosphate bridge but found {len(unique_pp_groups)}"
    
    # Check for a glycerol attachment indicator.
    # We require a branch with an 'OCC' fragment as a rough sign that the phosphate is linked to glycerol.
    glycerol_pattern = Chem.MolFromSmarts("OP(OCC)")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "Glycerol attachment not detected (missing typical OP(OCC) fragment)"
    
    # Count the acyl ester groups.
    # Acyl groups are represented by an ester pattern: OC(=O)[#6].
    acyl_pattern = Chem.MolFromSmarts("OC(=O)[#6]")
    acyl_matches = mol.GetSubstructMatches(acyl_pattern)
    if len(acyl_matches) < 2:
        return False, f"Found {len(acyl_matches)} acyl ester group(s); need at least 2"
    
    return True, "Contains cytidine headgroup with unique diphosphate bridge, glycerol attachment, and at least 2 acyl ester groups"

# If run as a script, test with one example.
if __name__ == "__main__":
    test_smiles = "P(OC[C@H]1O[C@@H](N2C=CC(=NC2=O)N)C(O)[C@H]1O)(OP(OC[C@H](OC(=O)CCCCCCC/C=C\\C/C=C\\CCCCC)COC(=O)CCCCCCCCCCCCCCC)(O)=O)(O)=O"
    result, reason = is_CDP_diacylglycerol(test_smiles)
    print(result, reason)