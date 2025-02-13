"""
Classifies: CHEBI:65102 3-hydroxy fatty acyl-CoA(4-)
"""
"""
Classifies: 3-hydroxy fatty acyl-CoA(4-)
Definition:
  "An acyl-CoA(4-) oxoanion arising from deprotonation of the phosphate and diphosphate OH groups 
   of any 3-hydroxy fatty acyl-CoA; major species at pH 7.3."
This classifier checks for:
  - A 3-hydroxy fatty acyl substructure as indicated by the presence of a hydroxylated carbon 
    two carbons away from a thioester (i.e., pattern "[C;H1](O)CC(=O)S").
  - A CoA moiety partial fragment (i.e., the "SCCNC(=O)CCNC(=O)" section commonly seen in CoA).
  - At least 4 negatively charged oxygens (O with formal charge -1) due to deprotonation.
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_3_hydroxy_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is a 3-hydroxy fatty acyl-CoA(4-) based on its SMILES string.
    
    The molecule should contain:
      - A 3-hydroxy fatty acyl chain (with the pattern [C;H1](O)CC(=O)S found in examples).
      - A CoA moiety as indicated by a fragment matching "SCCNC(=O)CCNC(=O)".
      - At least 4 deprotonated phosphate oxygens (O with a formal charge of -1).

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if molecule is a 3-hydroxy fatty acyl-CoA(4-), False otherwise.
        str: Explanation/reason for the classification.
    """
    
    # Parse the SMILES string to a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for 3-hydroxy fatty acyl part.
    # This pattern looks for a carbon with an -OH, followed by CH2 and a C(=O)S group.
    hydroxy_acyl_pattern = Chem.MolFromSmarts("[C;H1](O)CC(=O)S")
    if hydroxy_acyl_pattern is None:
        return False, "Failed to generate acyl substructure pattern"
    
    if not mol.HasSubstructMatch(hydroxy_acyl_pattern):
        return False, "No 3-hydroxy fatty acyl substructure pattern ([C](O)CC(=O)S) found"

    # Define SMARTS pattern for a CoA moiety fragment.
    # We require the fragment "SCCNC(=O)CCNC(=O)" which is common in all provided examples.
    coa_pattern = Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)")
    if coa_pattern is None:
        return False, "Failed to generate CoA moiety pattern"
    
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No CoA moiety fragment (SCCNC(=O)CCNC(=O)) detected"

    # Count the number of oxygens carrying a -1 formal charge.
    neg_oxygens = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 8 and atom.GetFormalCharge() == -1:
            neg_oxygens += 1
            
    if neg_oxygens < 4:
        return False, f"Insufficient deprotonated phosphate oxygens: found {neg_oxygens}, need at least 4"

    # Optional additional checks can include molecular weight, total atom counts, or rotatable bonds.
    # For this classification, our key substructure matches and negative charge count suffice.
    
    return True, "Molecule contains the 3-hydroxy fatty acyl chain and CoA moiety with at least 4 deprotonated oxygens (CoA(4-))"

# Example usage (uncomment for testing):
# test_smiles = "CCCCCCCCCCCCCC[C@@H](O)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OC..."
# print(is_3_hydroxy_fatty_acyl_CoA_4__(test_smiles))