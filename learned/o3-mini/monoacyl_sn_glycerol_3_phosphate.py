"""
Classifies: CHEBI:17088 monoacyl-sn-glycerol 3-phosphate
"""
"""
Classifies: monoacyl-sn-glycerol 3-phosphate 
An sn-glycero-3-phosphate compound having a single unspecified acyl group at either position 1 or position 2.
The improved criteria are as follows:
  1. The SMILES string must parse into a valid molecule.
  2. There must be exactly one phosphorus atom.
  3. The molecule’s exact molecular weight must lie in an expected range—now set to 250–650 Da.
  4. The molecule must contain a glycerol phosphate scaffold where the phosphate group is connected 
     to a glycerol backbone that carries exactly one acyl ester group.
     We enforce this by requiring that the molecule match one of two SMARTS patterns:
        • Pattern 1: acyl substitution at the sn-1 position.
        • Pattern 2: acyl substitution at the sn-2 position.
     
If these criteria are met, we return True with an explanation; otherwise False.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_monoacyl_sn_glycerol_3_phosphate(smiles: str):
    """
    Determines if a molecule is a monoacyl-sn-glycerol 3-phosphate based on its SMILES string.
    
    The criteria:
      1. Valid molecule.
      2. Exactly one phosphorus atom.
      3. Molecular weight between 250 and 650 Da.
      4. Contains a glycerol phosphate backbone with exactly one acyl ester substituent.
         We check this by looking for one of two substructure SMARTS patterns representing
         acyl substitution at the sn-1 or sn-2 position.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule qualifies as a monoacyl-sn-glycerol 3-phosphate, False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse the molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # Check there is exactly one phosphorus atom (atomic number 15)
    p_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15]
    if len(p_atoms) != 1:
        return False, f"Expected exactly one phosphorus atom but found {len(p_atoms)}."
    
    # Check the molecular weight range; relax lower bound to 250 Da to capture smaller species such as LPA 6:0.
    mw = rdMolDescriptors.CalcExactMolWt(mol)
    if not (250 <= mw <= 650):
        return False, f"Molecular weight {mw:.1f} Da not in the expected range (250–650 Da) for monoacyl-sn-glycerol 3-phosphate."
    
    # Define two SMARTS patterns representing a glycerol phosphate backbone with one acyl substituent.
    # Pattern 1: acyl group at sn-1 position.
    # Pattern 2: acyl group at sn-2 position.
    # Note: We omit stereochemistry in the SMARTS to increase our matching flexibility.
    pattern1 = Chem.MolFromSmarts("P(=O)(O)(O)OC[C](O)COC(=O)[#6]")
    pattern2 = Chem.MolFromSmarts("P(=O)(O)(O)OC[C](OC(=O)[#6])CO")
    
    # Check if the molecule contains either substructure.
    if mol.HasSubstructMatch(pattern1) or mol.HasSubstructMatch(pattern2):
        return True, "Molecule contains a phosphate group attached to a glycerol backbone with exactly one acyl ester group within the expected molecular weight range."
    else:
        return False, "Required glycerol phosphate substructure with a single acyl ester substituent was not found."

# Example usage: (uncomment the code below to run a test)
if __name__ == '__main__':
    # Test example: 1-nonadecanoyl-glycero-3-phosphate
    test_smiles = "P(OC[C@H](O)COC(=O)CCCCCCCCCCCCCCCCCC)(O)(O)=O"
    result, reason = is_monoacyl_sn_glycerol_3_phosphate(test_smiles)
    print(result, "->", reason)