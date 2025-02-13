"""
Classifies: CHEBI:17088 monoacyl-sn-glycerol 3-phosphate
"""
"""
Classifies: monoacyl-sn-glycerol 3-phosphate
An sn-glycero-3-phosphate compound having a single unspecified acyl group at either position 1 or position 2.
Improved criteria:
  1. The SMILES string must parse into a valid molecule.
  2. There must be exactly one phosphorus atom.
  3. The molecule’s exact molecular weight must lie in the expected range of 250–650 Da.
  4. The molecule must contain exactly one acyl ester group (the “OC(=O)” motif) – this enforces that there is only one acyl substituent.
  5. The molecule must have a glycerol–phosphate backbone.
     We check for a glycerol–phosphate “core” by matching a simplified SMARTS pattern.
     
If these criteria are met, we return True with an explanation; otherwise False.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_monoacyl_sn_glycerol_3_phosphate(smiles: str):
    """
    Determines if a molecule is a monoacyl-sn-glycerol 3-phosphate based on its SMILES string.

    The improved criteria:
      1. It must be a valid molecule.
      2. There must be exactly one phosphorus atom.
      3. Its exact molecular weight must lie in the range 250–650 Da.
      4. The molecule must contain exactly one acyl ester group – identified via the substructure "OC(=O)".
      5. The molecule must contain a glycerol–phosphate scaffold.
         We approximate this by searching for the substructure defined by:
           "P(=O)(O)(O)OC[C](O)CO"
         This pattern ignores stereochemistry and matches a phosphate attached to a glycerol backbone
         that still retains the second hydroxyl (the site that is not acylated).
         
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule qualifies as a monoacyl-sn-glycerol 3-phosphate, False otherwise.
        str: Explanation for the classification decision.
    """
    # 1. Parse the SMILES string into a molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # 2. Check that there is exactly one phosphorus atom.
    p_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15]
    if len(p_atoms) != 1:
        return False, f"Expected exactly one phosphorus atom but found {len(p_atoms)}."
    
    # 3. Check that the molecular weight is in the expected range.
    mw = rdMolDescriptors.CalcExactMolWt(mol)
    if not (250 <= mw <= 650):
        return False, f"Molecular weight {mw:.1f} Da not in the expected range (250–650 Da) for monoacyl-sn-glycerol 3-phosphate."
    
    # 4. Count the number of acyl ester groups.
    # We define an acyl ester group as the substructure "OC(=O)" (an oxygen attached to a carbonyl carbon).
    ester_pattern = Chem.MolFromSmarts("OC(=O)")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 1:
        return False, f"Expected exactly one acyl ester group (OC(=O)) but found {len(ester_matches)}."
    
    # 5. Check that the molecule has a glycerol-phosphate scaffold.
    # We define a simplified core as a phosphate (P(=O)(O)(O)) connected to a two-carbon chain ending in a hydroxyl,
    # which mirrors the glycerol backbone that still retains one free hydroxyl.
    # This pattern should be found when we remove the acyl chain from the glycerol.
    # For example, in "P(OC[C@H](O)COC(=O)...)(O)(O)=O" removing the acyl gives "P(OC[C@H](O)CO)..."
    glyco_pattern = Chem.MolFromSmarts("P(=O)(O)(O)OC[C](O)CO")
    if not mol.HasSubstructMatch(glyco_pattern):
        return False, "Glycerol–phosphate scaffold not found."
    
    # If all criteria are met, we assume the molecule is a monoacyl-sn-glycerol 3-phosphate.
    return True, ("Molecule contains one phosphate atom, has a molecular weight between 250–650 Da, " +
                  "contains exactly one acyl ester group (OC(=O)), and exhibits a glycerol–phosphate backbone.")


# Example usage (uncomment to test):
if __name__ == '__main__':
    test_smiles = "P(OC[C@H](O)COC(=O)CCCCCCCCCCCCCCCCCC)(O)(O)=O"  # 1-nonadecanoyl-glycero-3-phosphate
    result, reason = is_monoacyl_sn_glycerol_3_phosphate(test_smiles)
    print(result, "->", reason)