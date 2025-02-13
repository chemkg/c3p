"""
Classifies: CHEBI:65102 3-hydroxy fatty acyl-CoA(4-)
"""
"""
Classifies: 3-hydroxy fatty acyl-CoA(4-)
Definition:
  "An acyl-CoA(4-) oxoanion arising from deprotonation of the phosphate and diphosphate OH groups 
   of any 3-hydroxy fatty acyl-CoA; major species at pH 7.3."

The classifier now checks for:
  - A 3-hydroxy fatty acyl chain. Here we require that the molecule contains a carbon with a hydroxyl group 
    (written as [C@H](O) or [C@@H](O)) that is immediately followed by a CH2 and then by a thioester group 
    (i.e. the fragment should be –[C;H1](O)CC(=O)S–, which correctly locates the hydroxyl on carbon-3).
  - A CoA moiety fragment as indicated by the SMARTS fragment "SCCNC(=O)CCNC(=O)".
  - At least 4 deprotonated oxygens – since sometimes RDKit does not register the formal charge, we also count 
    occurrences of "[O-]" in the SMILES string.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_3_hydroxy_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is a 3-hydroxy fatty acyl-CoA(4-) based on its SMILES string.
    
    The molecule is required to contain:
      - A 3-hydroxy fatty acyl chain. We require a substructure where a carbon bearing one hydrogen 
        and an –OH ([C;H1]([OX2])) is immediately followed by a CH2 and then a thioester group (CC(=O)S).
      - A CoA moiety fragment ("SCCNC(=O)CCNC(=O)").
      - At least 4 deprotonated oxygens in the phosphate/diphosphate parts, determined either from the molecule’s
        formal charges or via matching the literal substring "[O-]" in the SMILES.
    
    Returns:
        (bool, str): A tuple of True and an explanation if the molecule meets our criteria,
                     or False and a reason if it does not.
    """
    
    # Parse the input SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a SMARTS pattern for the 3-hydroxy fatty acyl chain.
    # We require the fragment: a carbon with exactly one hydrogen and an attached hydroxyl, directly followed by a CH2 and then the thioester group.
    # This should match fragments like "[C@@H](O)CC(=O)S"
    acyl_pattern = Chem.MolFromSmarts("[C;H1]([OX2])CC(=O)S")
    if acyl_pattern is None:
        return False, "Failed to generate acyl substructure SMARTS pattern"
    
    if not mol.HasSubstructMatch(acyl_pattern):
        return False, "No clear 3-hydroxy fatty acyl substructure ([C;H1](O)CC(=O)S) found"

    # Define a SMARTS for the CoA moiety fragment.
    coa_pattern = Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)")
    if coa_pattern is None:
        return False, "Failed to generate CoA moiety SMARTS pattern"
    
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No CoA moiety fragment (SCCNC(=O)CCNC(=O)) detected"

    # Count the negatively charged oxygens by iterating over atoms.
    neg_oxygens_rdkit = sum(1 for atom in mol.GetAtoms() 
                             if atom.GetAtomicNum() == 8 and atom.GetFormalCharge() == -1)
    
    # Since RDKit may not always register the formal charges, we also count literal occurrences of "[O-]" in the SMILES string.
    neg_oxygens_smiles = smiles.count("[O-]")
    
    # Use the maximum of the two counts.
    neg_oxygens = max(neg_oxygens_rdkit, neg_oxygens_smiles)
    
    if neg_oxygens < 4:
        return False, f"Insufficient deprotonated phosphate oxygens: found {neg_oxygens}, need at least 4"
    
    # All requirements are met.
    return True, "Molecule contains the 3-hydroxy fatty acyl chain and CoA moiety with at least 4 deprotonated oxygens (CoA(4-))"

# Example usage (for testing):
# test_smiles = "CC\\C=C/C\\C=C/C\\C=C/C\\C=C/C\\C=C/C\\C=C/CC[C@H](O)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)n1cnc2c(N)ncnc12"
# print(is_3_hydroxy_fatty_acyl_CoA_4__(test_smiles))