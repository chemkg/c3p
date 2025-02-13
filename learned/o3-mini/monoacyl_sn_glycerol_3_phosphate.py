"""
Classifies: CHEBI:17088 monoacyl-sn-glycerol 3-phosphate
"""
"""
Classifies: monoacyl-sn-glycerol 3-phosphate
An sn-glycero-3-phosphate compound having a single unspecified acyl group at either position 1 or position 2.
"""
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

def is_monoacyl_sn_glycerol_3_phosphate(smiles: str):
    """
    Determines if a molecule is a monoacyl-sn-glycerol 3-phosphate based on its SMILES string.
    The improved criteria are as follows:
      1. The SMILES string must parse into a valid molecule.
      2. There must be exactly one phosphorus atom (to indicate a single phosphate group).
      3. The molecule must contain a glycerol backbone fragment.
         We use a simple heuristic SMARTS pattern "OCC(O)CO" (ignoring stereochemistry) for that.
      4. There must be exactly one acyl ester group.
         We detect this using the SMARTS pattern "OC(=O)[#6]" (an oxygen bonded to a carbonyl and then a carbon).
      5. The molecule’s exact molecular weight must lie in an expected range (roughly 300–650 Da)
         for monoacyl-sn-glycerol 3-phosphates.
         
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is considered a monoacyl-sn-glycerol 3-phosphate, False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse the molecule from the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # Check that there is exactly one phosphorus atom (for the phosphate group)
    p_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15]
    if len(p_atoms) != 1:
        return False, f"Expected exactly one phosphorus atom (phosphate group) but found {len(p_atoms)}."
    
    # Look for a glycerol backbone fragment.
    # We use a simple SMARTS pattern for glycerol: a three-carbon chain with two hydroxyl groups.
    glycerol_pattern = Chem.MolFromSmarts("OCC(O)CO")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "Glycerol backbone pattern (OCC(O)CO) not found in molecule."
    
    # Count acyl ester groups using a SMARTS pattern
    # "OC(=O)[#6]" should match an ester linkage of an acyl chain.
    acyl_ester_pattern = Chem.MolFromSmarts("OC(=O)[#6]")
    acyl_matches = mol.GetSubstructMatches(acyl_ester_pattern)
    if len(acyl_matches) != 1:
        return False, f"Expected exactly one acyl ester group but found {len(acyl_matches)}."
    
    # Check that the molecular weight is within a typical range (300–650 Da)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if not (300 <= mol_wt <= 650):
        return False, f"Molecular weight {mol_wt:.1f} Da not in the expected range (300–650 Da) for monoacyl-sn-glycerol 3-phosphate."
    
    # If all criteria are met, classify as monoacyl-sn-glycerol 3-phosphate.
    return True, "Molecule contains one phosphate group, glycerol backbone, and exactly one acyl ester group within the expected molecular weight range."

# Example usage:
if __name__ == '__main__':
    # Test example: 1-nonadecanoyl-glycero-3-phosphate
    test_smiles = "P(OC[C@H](O)COC(=O)CCCCCCCCCCCCCCCCCC)(O)(O)=O"
    result, reason = is_monoacyl_sn_glycerol_3_phosphate(test_smiles)
    print(result, "->", reason)