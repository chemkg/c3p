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
    It verifies that:
      1. The SMILES parses and that the molecule contains exactly one phosphorus atom,
         typical for a phosphate group.
      2. There is a glycerol backbone fragment, here heuristically defined as a 3-carbon chain
         with at least two hydroxyl groups (SMARTS: OCC(O)CO).
      3. The molecule contains exactly one acyl ester group (an oxygen linked to a carbonyl group)
         representing the single acyl chain. We search with the pattern "OC(=O)[#6]".
         
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is considered a monoacyl-sn-glycerol 3-phosphate, False otherwise.
        str: Explanation of the classification decision.
    """
    # Parse the molecule from the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Check for exactly one phosphorus atom (for the phosphate group)
    p_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15]
    if len(p_atoms) != 1:
        return False, f"Expected exactly one phosphorus atom (phosphate group) but found {len(p_atoms)}."
    
    # Look for the glycerol backbone.
    # Here we use a simple heuristic SMARTS pattern for glycerol: a three-carbon chain with two -OH groups
    glycerol_pattern = Chem.MolFromSmarts("OCC(O)CO")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "Glycerol backbone pattern (OCC(O)CO) not found in molecule."
    
    # Look for acyl ester group(s)
    # The pattern "OC(=O)[#6]" is used to find an ester linkage where an oxygen is linked to a carbonyl carbon.
    acyl_ester_pattern = Chem.MolFromSmarts("OC(=O)[#6]")
    acyl_matches = mol.GetSubstructMatches(acyl_ester_pattern)
    if len(acyl_matches) != 1:
        return False, f"Expected exactly one acyl ester group but found {len(acyl_matches)}."
        
    # Optionally, one may include additional checks (molecular weight, rotatable bonds, etc.)
    # For example, we can check that the molecule weight is in a reasonable range.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 250:  # approximate lower limit for such lipids
        return False, f"Molecular weight {mol_wt:.1f} Da is too low to be a monoacyl-sn-glycerol 3-phosphate."
    
    return True, "Molecule contains a phosphate group, glycerol backbone, and exactly one acyl ester group."

# Example usage (for testing purposes):
if __name__ == '__main__':
    test_smiles = "P(OC[C@H](O)COC(=O)CCCCCCCCCCCCCCCCCC)(O)(O)=O"  # 1-nonadecanoyl-glycero-3-phosphate
    result, reason = is_monoacyl_sn_glycerol_3_phosphate(test_smiles)
    print(result, "->", reason)